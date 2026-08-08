"""Tests for reroot_newick.py."""
import os
import subprocess
import sys
import tempfile
import unittest

try:
    from ete3 import Tree
except ImportError:  # pragma: no cover - depends which env the suite runs in
    Tree = None

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reroot_newick.py")

# Root is a polytomy with five children, which is the shape the real sampled
# trees have (HA/H7's root has five). set_outgroup must collapse it into a
# bifurcation without losing or duplicating a leaf.
POLYTOMY_ROOT = "((a1,a2)nodeA,b,c,(d1,d2)nodeD,e);"

# Target sits inside a polytomy rather than at the root.
NESTED_POLYTOMY = "(((p,q,r,s)nodeP,t)nodeI,u);"


@unittest.skipIf(Tree is None,
                 "ete3 not installed; run this module under envs/ete3.yaml")
class RerootNewickTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def run_script(self, newick, root):
        in_path = os.path.join(self.tmp.name, "in.nh")
        out_path = os.path.join(self.tmp.name, "out.nh")
        with open(in_path, "w") as handle:
            handle.write(newick)
        result = subprocess.run(
            [sys.executable, SCRIPT, "--input", in_path,
             "--output", out_path, "--root", root],
            capture_output=True, text=True,
        )
        return result, out_path

    def test_reroots_at_leaf_under_root_polytomy(self):
        result, out_path = self.run_script(POLYTOMY_ROOT, "c")
        self.assertEqual(result.returncode, 0, result.stderr)
        tree = Tree(out_path, format=1)
        # The outgroup leaf hangs directly off the new root, and the new root is
        # a bifurcation even though the original root had five children.
        self.assertEqual(len(tree.children), 2)
        self.assertIn("c", [child.name for child in tree.children])

    def test_preserves_leaf_set(self):
        for newick, target in ((POLYTOMY_ROOT, "c"), (NESTED_POLYTOMY, "q")):
            with self.subTest(target=target):
                result, out_path = self.run_script(newick, target)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(
                    sorted(Tree(out_path, format=1).get_leaf_names()),
                    sorted(Tree(newick, format=1).get_leaf_names()),
                )

    def test_preserves_unrooted_topology(self):
        """Rerooting must not change which splits the tree asserts."""
        result, out_path = self.run_script(POLYTOMY_ROOT, "c")
        self.assertEqual(result.returncode, 0, result.stderr)
        comparison = Tree(POLYTOMY_ROOT, format=1).compare(
            Tree(out_path, format=1), unrooted=True
        )
        self.assertEqual(comparison["rf"], 0)

    def test_target_inside_nested_polytomy(self):
        result, out_path = self.run_script(NESTED_POLYTOMY, "q")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("q", [child.name for child in Tree(out_path, format=1).children])

    def test_idempotent(self):
        """Rerooting an already-rerooted tree at the same leaf is a no-op."""
        result, first = self.run_script(POLYTOMY_ROOT, "c")
        self.assertEqual(result.returncode, 0, result.stderr)
        with open(first) as handle:
            once = handle.read()
        result, second = self.run_script(once, "c")
        self.assertEqual(result.returncode, 0, result.stderr)
        with open(second) as handle:
            self.assertEqual(handle.read(), once)

    def test_missing_target_exits_nonzero(self):
        result, out_path = self.run_script(POLYTOMY_ROOT, "not_in_tree")
        self.assertEqual(result.returncode, 1)
        self.assertIn("not_in_tree", result.stderr)
        self.assertFalse(os.path.exists(out_path))

    def test_internal_node_target_exits_nonzero(self):
        result, out_path = self.run_script(POLYTOMY_ROOT, "nodeA")
        self.assertEqual(result.returncode, 1)
        self.assertIn("internal node", result.stderr)
        self.assertFalse(os.path.exists(out_path))


if __name__ == "__main__":
    unittest.main()
