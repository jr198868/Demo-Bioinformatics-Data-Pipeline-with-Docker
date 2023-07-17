# Docker Image for Fastq Aata Process

This Docker image provides an environment for performing a bioinformatics data pipeline. It includes all the necessary dependencies and tools.


## Prerequisites

- Supported Operating System: Docker is compatible with a variety of operating systems, including Windows, macOS, and Linux. Make sure your system meets the requirements for running Docker. You can check the official Docker documentation for system requirements specific to your operating system.

## Installation Steps

Follow these steps to download and install Docker:

1. **Download Docker**: Visit the official Docker website (https://www.docker.com/) and navigate to the download page. Choose the appropriate version of Docker for your operating system and click on the download link to start the download process.

2. **Install Docker**:
   - **Windows**: Double-click on the downloaded installer file and follow the on-screen instructions to run the installation wizard. Accept the license agreement, choose the installation directory, and select any additional components you want to install. Once the installation is complete, Docker should be ready to use.
   - **macOS**: Double-click on the downloaded DMG file to open it. Drag and drop the Docker application into the Applications folder. Launch Docker from the Applications folder and follow any prompts to authorize the installation. Docker should be installed and ready for use.
   - **Linux**: Installation steps may vary depending on the Linux distribution you are using. Refer to the Docker documentation for detailed installation instructions specific to your distribution.

3. **Verify Installation**: After the installation is complete, open a terminal or command prompt and run the following command to verify that Docker is installed correctly:

```shell
docker --version
```

## Build a Docker image

To build the Docker image, navigate to the directory containing the Dockerfile and run the following command:

```docker
docker build -t bioinformatics_demo .
```

## Run a Docker container

To run a Docker container, you can use the following command.

Here's an example of how you can use the "-v" flag to mount a volume when running a Docker container:

```docker
docker run -it -d -v /path/to/local:/app bioinformatics_demo
```

This will execute the analysis pipeline inside the Docker container. Make sure to mount any necessary input files or directories to the container using the `-v` option.


## Copy files from the container to local machine folder

1. Check the running Docker containers on your system
```docker
docker ps
```

2. Find the 'CONTAINER ID' and run the following commands to copy analytical results to local machine folder

```docker
docker cp <CONTAINER ID>:/output/report.csv <local_path>
docker cp <CONTAINER ID>:/output/trim_read1_read2_alig.bam <local_path>
docker cp <CONTAINER ID>:/output/phred_quality_MAPQ_curve_read1.png <local_path>
docker cp <CONTAINER ID>:/output/phred_quality_MAPQ_curve_read2.png <local_path>
```

## Memory Allocation
To check the real-time resource usage statistics of running Docker containers, you can use the following command. It provides information about the CPU, memory, network I/O, and block I/O usage of each container.

```docker
docker stats
```


## Check execution time to run a Docker image

The execution time of the analysis pipeline depends on various factors, including the size of the input data and the complexity of the analysis. It is difficult to provide an exact estimate without specific details.

To determine how much time it takes to run a Docker image, you can use the time command in your shell.

```docker
time docker run <image_name>
```