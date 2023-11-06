# Template Repo

Template repository for data science projects.

## Introduction

The aim of this notebook is to play the role of a template for data science projects.

We are giving ourselves rules in order to improve code quality and generate a sharable and easy debuggable code.

## The workflow in a nutshell

Here we present the git workflow.

<p align="center">
 <img src="https://expressus.io/uploads/beautiful-gitflow-workflow-diagram.png">
</p>

Typically, what one does anytime wants to introduce an edit in a code is the following.

* Create a branch
* Add commits
* Open a pull request
* Discuss and review your code
* Implement a test pipeline
* If passed, merge.

In the previous session, we were sinners: we committed our code changes to our project on master branch.

In fact, generally the best approach to work with git is via branches and pull requests.

Of course, you can commit your code to master (Main / Production) branch but for avoid any mistakes and issues it is highly recommended to not do so.
You must create branch for any _task_ or _feature_ and then start working on that branch.

### Branches

A _branch_ is like a complete copy of your project but isolated.

#### Create a branch

When you create a branch in your project, you are creating an environment where you can try out new ideas without affecting the original code.
Indeed, changes you make on a branch do not affect `master` branch or any other branch, so you are free to experiment and commit changes.
You can act safely in the knowledge that your branch will not be merged until it is ready to be reviewed by someone you are collaborating with.
After that if you implemented a test pipeline, you can be reassured that code on other branches will be merged on master only if it passes all the tests.

__Note__: In most of real life scenarios in enterprises, you cannot directly push any commit toward master branch.
All changes must come with pull request.
In the most correct interpretation of such workflow, you actually have a second branch which is an exact copy of the master one, called `develop` or `staging` where you merge your feature branches on.
If pipelines pass on such branch than you merge this one on the master branch, that only accepts pull requests from this development branch.

Let's explore in a bit more detail some concept and related command.

#### Branching

Branching is a core concept in Git, and the entire GitHub flow is based upon it.
There is only one rule: anything in the master branch is always deployable.
Hence, keep the master branch clean and working.

Because of this, it is extremely important that your new branch is created off of master when working on a feature or a fix.
Your branch name should be descriptive (e.g., `refactor-authentication`, `user-content-cache-key`, `make-retina-avatars`), so that others can see what is being worked on.
Your commits should be small, frequent and they should have clearly explicative messages.

Usually you follow a convention for branch names and commit messages.
For example assigning to each branch, commit the corresponding Jira ticket code (or any other tracking system, like issues numbers in GitHub).

A commit must be self-contained and implement a single change in the codebase.
Do not hesitate to split changes into multiple commits if needed.

For instance, if you start a feature issue and end up writing a utility, fixing a typo in the doc, updating a test dependency and then implementing your feature, this would be 4 commits:

* `feat: Add the fooToBar utility, ABC-1`
* `docs: Fix typo in the Squirrel README, ABC-1`
* `test: Update the Foo dependency`
* `feat: Implement new Foo, ABC-1`

If a commit has nothing to do with your original ticket, it should ideally be done later in another branch.
Sometimes, that is not possible (for instance if someone force-pushed a non-linted file and your linter will not let you push, or if a test has a broken dependency).
If that happens, the commit message may not need the ticket suffix as it is a general maintenance commit.

##### Prefixes

The prefix qualifies the content of the commit , eg. "feat" to describe a new feature, "fix" for a bug fix, etc.
Each prefix is linked with a versioning strategy used by semantic-release to increment (or not) the project version.
For instance, “test” commits do not trigger the release of a new version.
Example of prefixes are

* feat: A new feature, ie. something that changes the behaviour of the code
* fix: A bug fix in the code (warning: fixes to the tooling should use chore)
* docs: Documentation only changes
* style: Changes that do not affect the meaning of the code, eg. white-space, formatting, missing semi-colons (warning: do not use this for CSS changes! CSS changes are feat)
* refactor: A code change that neither fixes a bug nor adds a feature, eg. moving a piece of existing code to a new util or renaming a class
* perf: A code change that improves performance
* test: Adding missing or correcting existing tests
* chore: Changes to the build process or auxiliary tools and libraries such as documentation generation

###### Example

* `feat(header): Turn the color of the header bar into blue`
* `style(auth): Harmonise whitespace in the auth module code`

#### Commit message relevance

The commit message has to precisely describe, in English, the change brought by the commit, as follows.
The first word is an action verb in imperative form.
It is capitalised. No punctuation should be used.
The message should be short.

```bash
prefix(scope): Message starting with action verb
```

_e.g._: `feat(header): Turn the color of the header bar into blue`

#### Git Commit Message Hook

There is a [python package](https://git-message-hook.readthedocs.io/en/latest/) that can be installed and help to force commit messages in the format described above.

The only step necessary is to install such a package by

```bash
pip install git-message-hook
```

This will automatically add a commit message hook in the .git/hook directory, make it executable by

```bash
chmod +x .git/hook/commit-msg
```

and that's it.

#### Branch naming

Branch type.
Following the GitHub branch model, each branch should be prefixed with a specific type.

##### Feature branch

Used for specific feature work or improvements. Generally branches from, and merges back into, develop, using pull requests.

`feature/xxx`

###### Bugfix branch

Typically used to fix develop.

_e.g._: `bugfix/xxx`

##### Hotfix branch

Used to quickly fix production/staging branches without interrupting changes in develop.

Changes made in a hotfix branch have to be merged back to develop once merged into production or staging.

_e.g._: `hotfix/xxx`

In order to respect such rules one may want to install a [pre-commit hook](https://pre-commit.com/), _i.e._ an automatic check that your commit respects all the rules above.

---

#### Try the following

See list of branch and the active one with an asterisk (*) on the left.

```bash
git branch
```

Create a new branch

```bash
git branch [branch name]
```

Switch to the newly created branch :

```bash
git checkout [branch name]
```

You can do that in one shot.
Create a branch and switch to it.

```bash
git checkout -b [Branch name]
```

#### Add commits

What you want to do is making changes in the new branch environment and committing your code!

You are basically creating a transparent history of your work in the newly created branch isolated from master branch so do whatever you want and commit your code!

After changes you need to push your branch to GitHub.

```bash
git push origin <Name of your branch>
```

#### Pull request

You are now ready to your first pull request.
Create a pull request to propose and collaborate on changes to a repository.
These changes are proposed in a _branch_, which ensures that the `master` branch only contains finished and approved work.

If you check github repository now, you can see new branch in Web GUI and you can easily Click on `Compare & Pull request`.

You can put see changes and compare them with master or any other branches.
You can comment and mention other team members to review your code, basically they can see all the modifications and all the changed files, make suggestions, etc.

[For more information have a look here](https://help.github.com/en/articles/creating-a-pull-request).

By pipelines or _actions_ in GitHub, you can configure continuous integration but for now we are just going to stick with default settings.

Once you have created a pull request, you can push commits from your topic branch to add them to your existing pull request.
These commits will appear in chronological order within your pull request and the changes will be visible in the "_Files changed_" tab.

Other contributors can review your proposed changes, add review comments, contribute to the pull request discussion, and even add commits to the pull request.

##### Merge

After you are happy with the proposed changes, you can _merge_ the pull request.
If you are working in a shared repository model, the proposed changes will be merged from the head branch to the base branch that was specified in the pull request.

##### Rebase

The last concept we want to explore is the _rebase_ of a branch onto another.

Imagine you have created a new branch from master and you are working on a feature.
Your colleague also is implementing a fix, so he created another branch.
He was quicker than you, so he merged to master and you are still working on your branch, that is not aligned to the master anymore (as your colleague merged his fix).
If you open a pull request, you will find that this would be impossible to merge.
How can you solve such an issue?

The figure below helps us.
By a rebase, you copy all the commits that have been made on the base branch after you created yours in the correct order.

<p align="center">
    <img src="https://docs.gitlab.com/ee/topics/git/img/git_rebase_v13_5.png">
</p>

After a rebase you can ask your pull request to be merged!

---

## Pre-commit hook

Before submitting your code to code review, a useful tool to identify issues are pre-commit hooks.
Indeed, by pointing these issues out before code review, this allows a code reviewer to focus on the architecture of a change while not wasting time with trivial style nitpicks.

Pre-commit scripts are usually in a hided directory of the repository, and you give them the permission to be executed on your machine.
Before you can run hooks, you need to have the pre-commit package manager installed.

Using pip:

```bash
pip install pre-commit
```

In a python project, add the following to your `requirements.txt` (or `requirements-dev.txt`):

`pre-commit==3.0.xxx`

For further details, have a look at [pre-commit.com](https://pre-commit.com/).
