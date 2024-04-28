<!-- Adapted from d-SEAMS: https://github.com/d-SEAMS/seams-core -->
# Contributing

We love pull requests from everyone. By participating in this project, you
agree to abide by the contributor covenant suggested [code of conduct].

[code of conduct]: https://github.com/flowy-code/flowy/blob/main/CODE_OF_CONDUCT.md

Do make sure to run tests and generally **BE PREPARED** to have your code vetted and
checked. Do **NOT** submit code you would not be prepared to defend or maintain,
unless you are fixing a bug.

Push to your fork and [submit a pull request][pr].

[pr]: https://github.com/flowy-code/flowy/compare

At this point you're waiting on us. We like to at least comment on pull requests
within four business days (and, typically, three business day). We may suggest
some changes or improvements or alternatives.

Some things that will increase the chance that your pull request is accepted:

- Write tests.
- Follow the commit-style below.

## Commit Style

A sample **good commit** is:

```diff
fileName: Thing I did
Some subHeading things

So this change was to do that thing I thought was good. Also there was this
other person who thought so too, so then I ate a sandwich and we got the code
done. I am writing this but really, honestly, two people did this.


Co-authored-by: Joel Doe <joel@iexistreally.com>
```

A good commit should have:

- The name of the file, or the topic or the subject you have changed or the
  namespace or the functionality you have added **something:**
- A line describing **something:**
- An _(optional)_ subheading with more details
- An _(optional)_ paragraph or essay on why the change was done and anything else you want to share with the devs.
- **Co-authorship** IS MANDATORY if applicable. Even if you just had a sandwich
  with the other person. It won't kill you to share, or to write that.

<!-- * Follow our [style guide][style]. -->

<!-- [style]: https://github.com/thoughtbot/guides/tree/master/style -->

### Commit template
If you are not already using a commit template, consider doing so. The Seldon repository includes a template under `.gitmessage`.

```sh
git config --global commit.template .gitmessage
```

You can omit the `--global` to only use the template when committing to Seldon.

# Automated Styles

<!-- Kanged from readcon: https://github.com/HaoZeke/readCon -->

A `pre-commit` job is setup on CI to enforce consistent styles, so it is best to
set it up locally as well (using [pipx](https://pypa.github.io/pipx) for
isolation):

```sh
# Run before commiting
pipx run pre-commit run --all-files
# Or install the git hook to enforce this
pipx run pre-commit install
```

# Changelog management

We use markdown for the changelog fragments, and they are stored under `docs/newsfragments`

- Use issues to number the snippets
- Build a final variant with: `towncrier build --version 0.6.3 --date "$(date -u +%Y-%m-%d)"`
- Supported categories are:
 + `security`
 + `removed`
 + `deprecated`
 + `added`
 + `changed`
 + `fixed`

Here are some sample use cases:
```sh
towncrier create -c "Fancy new feature but without an issue attached" +new_feat.added.md
towncrier create -c "Require C++17 only" 1.changed.md
```

The generated markdown files can be modified later as well.
