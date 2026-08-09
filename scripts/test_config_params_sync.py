"""Guard the one place a config value can go stale without anyone noticing.

Most rules hand their script a resolved value (`--root '{params.new_root}'`),
so what Snakemake hashes for the params rerun-trigger is literally what the
script uses, and the two cannot drift. download_all_references is the
exception: it passes `--config config.yaml` and the script re-reads the file
itself, so the rule's `params:` block is only Snakemake's *belief* about which
values matter. Add a key to the script without adding it to the rule and
editing that key silently rebuilds nothing -- the exact failure mode removing
config.yaml-as-an-input was meant to eliminate, relocated somewhere harder to
see.

This asserts the belief matches the code.
"""
import os
import re
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "scripts", "download_ref_seq.py")
SNAKEFILE = os.path.join(REPO, "Snakefile")


def config_keys_read_by(path):
    """Top-level config keys a script reads, from `config["key"]` occurrences."""
    with open(path) as handle:
        return set(re.findall(r'config\[["\'](\w+)["\']\]', handle.read()))


def rule_params_block(name):
    """The text of one rule's `params:` block, up to the next top-level field."""
    with open(SNAKEFILE) as handle:
        text = handle.read()
    start = text.index(f"rule {name}:")
    end = text.index("\nrule ", start + 1)
    body = text[start:end]
    params_at = body.index("\n    params:")
    rest = body[params_at + 1:]
    # Fields are indented four spaces; params entries are indented eight.
    following = re.search(r"\n    (?!    )\w+:", rest)
    return rest[:following.start()] if following else rest


class DownloadReferencesParamsSyncTestCase(unittest.TestCase):
    def test_every_config_key_the_script_reads_is_declared_as_a_param(self):
        params = rule_params_block("download_all_references")
        missing = sorted(k for k in config_keys_read_by(SCRIPT)
                         if f"config[\"{k}\"]" not in params
                         and f"config['{k}']" not in params)
        self.assertEqual(
            missing, [],
            "download_ref_seq.py reads these config keys, but "
            "download_all_references does not declare them in params:, so "
            "editing one would not re-trigger the download: " + str(missing),
        )

    def test_the_script_actually_reads_some_config_keys(self):
        """Guards the guard: a regex that matches nothing would pass silently."""
        self.assertTrue(config_keys_read_by(SCRIPT))

    def test_config_yaml_is_not_an_input_of_the_rule(self):
        """If it becomes one again, this whole test stops being necessary."""
        with open(SNAKEFILE) as handle:
            text = handle.read()
        start = text.index("rule download_all_references:")
        body = text[start:text.index("\nrule ", start + 1)]
        inputs = body[body.index("\n    input:"):body.index("\n    params:")]
        self.assertNotIn("config.yaml", inputs)


if __name__ == "__main__":
    unittest.main()
