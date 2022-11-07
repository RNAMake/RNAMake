# RNAMake Development Guide

This document is a guide for users wishing to make code contributions to `RNAMake`.

## Initial Steps

To get started, you first need a GitHub account and `git` installed on your local machine. If you need help installing `git`, follow along with these guides:
* [How to install git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
* You may have to follow along with these authentication guides as well:
  * [Generating an SSH Key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
  * [Creating a personal access token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)

You also need the following tools installed on your computer:
* a C++ compiler (GCC, Clang, etc.) and CMake
* an IDE ([CLion](https://www.jetbrains.com/clion/promo/) is free with a student license)
* Docker (You can follow this [Docker installation guide](https://docs.docker.com/engine/install/) if necessary)

**_NOTE:_** You only need the compiler and CMake if you want to develop locally. Alternatively, you can develop against the Docker container, which we'll cover later.

Once you have everything you need, clone the project with the following commands:

```
git clone git@github.com:RNAMake/RNAMake.git
```

You should now have a directory called `RNAMake` that contains all the RNAMake code.

## Development Workflow

The `RNAMake` project follows a trunk-based development model. This means all changes start by making a new branch from the `development` branch, making small changes, and then creating a pull request. Let's discuss each step in detail:

### Creating a feature branch

All development will happen against the `development` branch. This means that before you write any code, you need to first checkout `development` and pull any changes that may have happened since you last pulled the branch. You can do this by running the following commands (note the `$>` represents your command prompt and is not part of the command to be run):

```
$> git checkout development
$> git pull
```

This will synchronize your local copy of the `development` branch with the remote `development` branch (the one that lives on GitHub). The last step to take before you start working on your changes is to create a new branch called a "feature branch." For this guide, we'll call the branch `my-new-feature`. You can optionally prepend your initials and a `/` to the branch name; for example, I might call my branch `ew/my-new-feature`. To create this new branch, first make sure you're still "on" the development branch, you can do this by running:
```
$> git branch
```
And you should see the word `development` either with an asterisk next to it or colored differently from the other branches. Once you're sure you're on the development branch, we can now create a feature branch based on `development`. This is accomplished with the following command:

```
$> git checkout -b my-new-feature
```

You can double check that you're on the `my-new-feature` branch by running `git branch` again and seeing `my-new-feature` indicated either with an asterisk or different color. Once you're on a feature branch, you're ready to start actually writing code!

### Developing with Docker

Now that you're ready to start writing code, you should know that `RNAMake` supports Docker development. There is a verbose guide for how to develop with Docker in the RNAMake project in this repository already:

[RNAMake Docker Development Guide](help_with_docker.md)

You're free to develop locally without Docker if you want but be advised, part of the CI pipeline (more on that later) checks that the Docker images can be built. If they cannot, your pull request will not be approved until you fix the issues.

Please also note that our Docker development process is new and is undergoing changes often. Please check Slack or the [RNAMake Docker guide](help_with_docker.md) often to see if there are updates to the process. As of right now, the process is not optimal and Erik is working on making it better. In the meantime, feel free to reach out to Erik with any questions.

### Commiting Your Work

Once you've made all the code changes you want to make, it's time to commit your work. This means you will make your changes an official part of your feature branch's history. You can think of commits as save points that you can come back to later if things go awry. In order to make a commit, you first have to add the files that you changed:

```
$> git add filename1 filename2
# OR
$> git add .
```

This sets up the files you updated to be "tracked." Now, you actually make the commit. Think of a good message that summarizes the work you're commiting and run the following command:

```
$> git commit -m "My commit message"
```

This will add your code changes to the commit history of your branch, meaning they will show up as changes when you make a pull request against the `development` branch.

As a shortcut, you can add and commit changes in one shot by passing the `-a` flag to `commit` like this:

```
$> git commit -am "My commit message"
```

However, for some reason, passing the `-a` flag will not add newly created files. So if you created a brand new file that was not part of the repository before, you'll most likely have to run `git add .` before commiting.

### Publishing your changes

Once you've commited all the work you're doing, it's time to publish your changes and create a pull request. You do this by running the following command:

```
$> git push --set-upstream origin my-feature-branch
```

(Note that everything after `push` only has to be done the first time you push a branch)

This command will publish your branch to the remote repository (that is, to GitHub). From there, you can make a pull request. A pull request is a request for the maintainer of the project to add your code changes to the project. You create a pull requeste (or just "PR") by going to the remote repository (https://github.com/RNAMake/RNAMake), clicking "Pull Requests" in the reposiotry menu, clicking the "New pull request" button, and selecting which branches you want to be part of the PR. Usually, the branch on the left should be `development` and the branch on the right will be your new branch, `my-new-feature` in our previous examples.

Often, if you go to the repository website shortly after pushing your branch, there will be a prompt for creating a pull request based off your new branch. You can click that "Create pull request" button too, it is the same process.

On the PR page, you can write about your changes before clicking "Create pull request" button in the bottom right. You do need to click that button to complete the process of making the PR though.

### Pull Request Checks

Every time there is a pull request made, GitHub automatically runs three checks: continuous integration build, a code quality check, and an approval from either Dr. Yesselman or Erik:

#### Continuous Integration (CI) with CircleCI

The CircleCI build process (which we just call "the build") is very important. The first thing it does is build a Docker image based off of your code changes, then runs the image. Finally, inside that Docker image, it runs all the unit tests in the project. If any of these steps returns an error code, the build has failed and you must fix whatever caused the fail to happen before the code will be merged.

**_The build MUST PASS before pull requests are merged!_**

If you build fails for some reason, you can go to [the CircleCI dashboard](https://app.circleci.com/pipelines/github/RNAMake) and click on your build to see the reason the build failed. The CircleCI free-tier only allows a few users to have access to the project however, so only a few lab members will be able to see the CircleCI build.

However!

Because we are using Docker development, you should not ever need to access CircleCI since building the Docker image and running tests against the Docker container will be the same no matter what environment they're ran on. So 99% of the time your build breaks, you can reproduce the problem by running:

```
$> docker build -t rnamake -f docker/Dockerfile .
$> docker run -d -i --name RNAMake rnamake
$> docker exec rnamake_dev ctest -V --test-dir /RNAMake/cmake/build
```
If any of those steps fail on your local computer, they will also fail in the build, so it is best to make sure these steps pass for you before making a pull request.

#### Code Quality Check with CodeClimate

CodeClimate is a tool that parses the changes you've made to your PR for common code-quality type problems. Things like creating a variable and then never using it, incorrect indentation, lines that are too long and more are all things that CodeClimate is configured to catch. CodeClimate is another tool which limits how many people can access it for free accounts. However, CodeClimate (or just CC) is configured to leave comments on pull requests, so you should never have to log into CC.

Sometimes (often), CodeClimate will complain about something that you don't really have a good way to fix. For example, CC will often complain that a file is too long, but this is the kind of thing that very often not worth fixing. In such cases, the PR reviewer can mark those items as `wontfix` and you shouldn't have to worry about fixing something liek that. However, you should fix any legitimate complaints CC has about your code. CC is one of the three PR checks so you can't merge your changes until its step passes, just like with CircleCI.

#### Code review

Every pull request is subject to "code review." This is basically where another developer looks over your changes and makes comments about them. The developer may ask you to change something, why you did something, or just compliment your work so far. Either way, the final check for a PR is that a reviewer must "approve" the pull request. Right now, the only people allowed to approve pull requests are Dr. Yesselman and Erik (Erik will only approve PRs that touch things like the Dockerfile or CI config file). You can request a review by clicking "Reviewers" in the right column of the PR page, or you can simply @ them and request a reivew (for example, you can leave a comment that says "@erik-whiting can you please review and approve this PR?" and Erik will be notified).

### Merging

Once you've made your changes and all the checks have passed, it's finally time to merge your work into the `development` branch, congratulations! Only Dr. Yesselman and Erik have PR merging privileges so you will have to ping one of them when your PR is ready to go.

This part is important: Once your PR is merged into `development`, your (and everyone else's) local `development` branch will be out of sync with the remote `development` branch. To sync your local with the remote `development` branch, you should checkout `development` and `pull` the changes:

```
$> git checkout development
$> git pull
```

This will resync the two branches and you'll be ready to start the process all over again.

## Do you have questions? Did you break something?

If you have questions or if something goes terribly wrong, please reach out to Erik either on Slack or via email (ewhiting4@huskers.unl.edu). In most cases, he will be able to answer any of your git, GitHub, Docker, and CircleCI questions. Don't be shy, that's what he's here for!
