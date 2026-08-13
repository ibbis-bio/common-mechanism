# Contributing to Commec

Thanks for your interest in contributing to Commec. Commec is maintained
by a team of IBBIS staff, but we welcome external contributions. If you
are outside IBBIS and planning a non-trivial change, please open an issue
first so we can discuss the approach before you invest time in a PR. Note
that Commec is [MIT-licensed](/LICENSE).

This doc describes our branching process and PR checks.

## Branching model

We use a version of [Gitflow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow)
with a main branch containing only commits that are ready to be part of tagged versions.

- **`main`**: should only contain stable, release-ready code (lets others clone
  the repo without pulling something that isn’t fully validated). Do not target
  `main` with feature work.
- **`develop`** — should be the starting point for new feature branches, i.e.
  approved + tested code that isn’t yet ready for release
- **`release_v#.#.#`**: branches created from `develop` to collect a release’s
  worth of commits together for functional testing before a merge to main

```
git checkout develop
git pull
git checkout -b my-feature      # descriptive name, e.g. fix-cmsearch-evalue
```

You can name your feature branches whatever you want. They should be squash merged
into `develop` (In the rare case where a feature branch contains important commit history, use an interactive rebase to clean the history + and rebase merge those commits onto develop)

## Before you open a PR

Run the same checks CI will run. Set up the environment once:

```
conda env create -f environment.yaml
conda activate commec-env
```

Then, from the repo root:

```
ruff check --select I --fix    # import sorting
ruff format                    # formatting
pytest -vv                     # test suite
```

CI (`Test Commec`, defined in `.github/workflows/pr-checks.yml`) runs the same
ruff checks and the full `pytest` suite against every PR to `develop` and
`main`. PRs must be green to merge. Running locally first saves a round-trip.

## Databases

Commec's reference databases (biorisk, control lists, best-match, low-concern)
are versioned and released **separately** from this package, and are distributed
via Cloudflare R2 rather than living in this repository. `commec setup` fetches
them. If your change affects database compatibility, update
`requires-databases.yaml` accordingly.

## Questions

Open an issue, or reach the maintainers through the channels listed in the
repository README.