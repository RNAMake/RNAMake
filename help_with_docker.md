# What's Docker?

According to the Docker docs:
> Docker is an open platform for developing, shipping, and running applications. Docker enables you to separate your applications from your infrastructure so you can deliver software quickly. With Docker, you can manage your infrastructure in the same ways you manage your applications. By taking advantage of Docker’s methodologies for shipping, testing, and deploying code quickly, you can significantly reduce the delay between writing code and running it in production.

Basically, with Docker, we can run a containerized system on our own computers. For example, if I want to run RNAMake in a Linux environment, but all I have is Windows, I can make a *Dockerfile* that contains a definition of how to build such a Linux system as well as the compile steps (see the [devDockerfile](docker/devDockerfile) in this repo) I would run on a real Linux machine. If I run this image, Docker creates a *container* based on that image. It's like running a lightweight virtual machine.

> Note: For some reason, people can get really upset if you refer to Docker as a "lightweight VM" so don't use that term outside of the lab.

## Basic Usage

In the context of this project, the first thign we'll do is create a development image. This is done by running the `build` command and passing `docker/devDockerfile` to it. In your terminal, `cd` into the RNAMake directory then run the following:

```sh
$> docker build -t rnamake_dev -f docker/devDockerfile .
```

This will build the development image. Basically, Docker will spin up a CentOS machine, install all the things listed in that file, and then build and compile the project. The build takes anywhere between 5 and 10 minutes. Please note, the `-t` flag stands for "tag" and is what you name the image. You can call it whatever you want.

### Using `docker-compose`

The `docker-compose` command line tool helps coordinate your docker-based applications. You use it by creating a `docker-compose` file that configures your images, volumes, and so on. In this project, once you've created the `rnamake_dev` container, you can run the following command in the terminal:

```sh
$> docker-compose up -d
```

This will start the docker container and will also mirror all your local files with the files inside of the container. This means that you can update code files locally in your IDE and expect them to show up in the container without having to rebuild the whole image.

To recompile, you can use the following command:

```sh
$> docker exec rnamake_dev ninja
```

To run the test suite, run:

```sh
$> docker exec rnamake_dev_ ctest
```

### Using Traditional `docker-run`

Once the image is done building, you then have to *run* the image before you can do anything with it:

```sh
$> docker run -d -i --name RNAMake_dev rnamake_dev
```

Note that the last argument is whatever you named the image in the `build` command. The `-d` flag tells Docker to run the image in "detach" mode. This means the image is run in a background process; this is necessary if you want to continue using your terminal after running the container. The `-i` command means "interactive" and it tells Docker to keep the image running after it's done building. This is equivalent to keeping the container "on" once it's done running. If you forget this flag, the container will start up and then immediately turn itself off, preventing you from using it.

Speaking of using it! Once the image is running, you can get into it a few different ways. If you have Docker desktop, you can just click on the container and click "Open in terminal." Alternatively, you can get into the container by running the following command:

```sh
$> docker exec -i -t RNAMake_dev bash
```

Note that the second to last argument is the name of the container from your `run` command. The above command tells Docker to execute the command "bash" against the container `RNAMake_dev` and keep the connection open. `bash` is the command line utility for most Linux distros, so this is similar to SSHing into a server.

### Did you get an error on the `build` command?

If you run into an error when running the build command that says something about the Docker daemon not running, you have a few options. If you have Docker desktop (which I recommend), just run it, that will start the daemon.

If you don't have Docker desktop, there are a variety of ways to start the daemon depending on what operating system you're on. The Docker license has changed in recent years so a lot of guides online are out of date. As such, you might have to try a few things. If you spend more than 20 minutes trying to get the daemon started though, reach out to Erik (ewhiting4@huskers.unl.edu) and he will help you get up and running.

If you got an error but it wasn't about the Daemon, it might be because you're on a mac. If so, you might have to run

```sh
$> docker build -t rnamake_dev --platform linux/amd64/v8 -f docker/devDockerfile
```

# To build a minimized image (i.e., just RNAMake functionality)

The following commands are similar to those listed above, except they use `Dockerfile` instead of `devDockerfile`. The difference between this Dockerfile and the development once, is that this one *does not* keep the source code files around. This means we can share this image with other people without worrying about them stealing our code.
```sh
$> docker build -t rnamake -f docker/Dockerfile .
$> docker run -d -i --name RNAMake rnamake
```

# Ubuntu Image

If you don't like the CentOS dev image, there's an Ubuntu one as well. You can build it with

```sh
$> docker build -t rnamake_ubuntu -f docker/ubuntuDockerfile
```

Note that this is a development image--the source code files are persisted in this one--so don't share this one with anyone outside of the lab.

# Still need help?

Get a hold of Erik (ewhiting4@huskers.unl.edu) if you're stuck and need assistance.
