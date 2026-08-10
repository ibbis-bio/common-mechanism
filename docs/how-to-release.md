# Commec Release Instructions.

What branching system do we use?
We use a version of Gitflow with a main branch containing only commits that are ready to be part of tagged versions. We don't squash commits into main. Feature branches are squash merged into Develop. When a release is ready, a release branch is made from Develop, final tests and changes are made, and it is merged into Main, and merged back into develop.

### Branches:

**main**: should only contain stable, release-ready code (lets others clone the repo without pulling something that isn’t fully validated)

**develop**: Consists of approved code that isn’t yet ready for release. This should be the starting point for new feature branches.

**feature branches**; you can name these whatever suits; they should be squashed + merged into develop (unless they have a ton of commits, in which case you may wish to instead use an interactive rebase to clean up the history + a merge commit).

**release_v#.#.#**: Branches created from develop to collect a release’s worth of commits together for functional testing before a merge to main.

>Note: Do not rebase develop, do not rebase main, this will break true CI/CD flow. As altering the history of these branches will affect how all features and releases are merged, and create conflict.

How do you create a new release of the commec package?

## 1. Create and test a release branch
The time has come to create a new commec release! Create a branch from develop, and open a draft PR against main. The branch should be named “release_v#.#.#”.

Do any last minute functional testing. Any commits made in response to testing should be ideally be squashed into a single commit.

Update the version in code (i.e. in pyproject.toml, meta.yaml, ensure the json functional test output version is also updated).

>Now is also a good time to ensure that the conda recipe yaml is as up to date as possible. Note you cannot supply the sha256 code, as this results in a chick-egg situation. The repo conda recipe is a reference, and not the source of truth. Dependencies are ultimately stored in the environment.yaml.

## 2. Merge release branch into main, and develop.
Merge the release branch into main,  bringing in the develop commits + version increment commit. The commit will look like this one: 4ebd5a2

Also immediately merge the release branch back into develop. This may result in conflicts based on how much further developement has occured in main since the release branch occured, but will be crucial to reduce conflicts in the future.

## 3. Create a new tagged version in Github
Go to https://github.com/ibbis-bio/common-mechanism/releases and click “Draft a new release”

Create a new tag. We follow semantic versioning, as follows:

* **Major release** (e.g. 0.x.x → 1.x.x): breaking changes, no promise of backwards-compatibility (e.g. complete changes to output format)
* **Minor release** (e.g. 0.1.x → 0.2.0): backwards-compatible changes that include new features
* **Patch** (e.g. 0.1.0 → 0.1.1): bugfixes, refactors, other small and backwards-compatible changes that don’t include new features

Write a description of the release, following this format (make sure to update the version):
```
Name: commec v2.0.0

Description:

## Features
* Add a human-readable description of the features. Feel free to include direct links to PRs.

## Bug fixes
* * Add a human-readable description of the bugfixes. Feel free to include direct links to PRs.

**Full Changelog**: https://github.com/ibbis-screening/common-mechanism/compare/v0.1.2...v0.2.0
```
## 4. Create a new version of the bioconda recipe
The easiest method for users of commec to install is with a conda package manager. To be able to continue using this method with the latest version, the Commec recipe requires updating.
### 4.1: Updating meta.yml
Go to the [bioconda recipe fork](https://github.com/ibbis-bio/bioconda-recipes), and "sync fork", discard any existing commits, it is easiest to work from the latest version of bioconda recipes.
Update the recipes/commec/meta.yaml file, lines 2 (version number) and 3 (sha256) of the meta.yml
The Version number should match the version number of the newly cut release. Format is X.X.X
The shasum256 should be generated using ```sha256sum filename``` (for linux) and pasted exactly. The shasum is for the tar.gz file that is created as part of the tag release. It can be downloaded, and sha256sum calculated.
Add/remove any updates to list of dependencies, (should be same as environment.yml in Commec source)

If for some reason the commec version is the same but another conda recipe release is occuring (Updating the bioconda recipe but not commec, due to a broken test, dependency, or updated sha256sum), simply increment the field:
```
build: 
   number: 0
```

### 4.2: Local Testing and building (OPTIONAL)
Test the build locally (recommended in a new conda environment)  using 
```
bioconda-utils lint –packages commec
bioconda-utils build –packages commec
```
Once the build passes all tests, we're ready to merge to bioconda.

>Note: This process can be surprisingly quite memory intensive, straining 16 GB, and the process may kill itself with little or no error logging in this case.

>Note: You could try `bioconda-utils --autobump` see the [bioconda-utils command line reference](https://bioconda.github.io/contributor/bioconda-utils-cli.html) which can manage the creation of the new pull request etc, and gets update time automated to 1-3 days.

### 4.3: Deployment
Create a Pull request for the forked bioconda-recipes repository. Easily done in browser with the `contribute` button on the forked repo landing page. (next to `sync fork` from earlier). Name the PR "Update to Commec vX.X.X". The PR will undergo automated deployment testing. 

If anything fails, you can simply commit changes to the forked repo, and it will redo automated testing.

Once testing is successful, comment 
```“@BiocondaBot please add label”```, and in 1-7 days time, someone will approve the PR.

Commec updated!
