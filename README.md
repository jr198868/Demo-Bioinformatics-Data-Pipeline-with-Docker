# Docker Image for XYZ Analysis

This Docker image provides an environment for performing XYZ analysis. It includes all the necessary dependencies and tools.

## Build a Docker image

To build the Docker image, navigate to the directory containing the Dockerfile and run the following command:

```docker
docker build -t bioinformatics_demo .
```

## Run a Docker container

To run a Docker container, use the following command:

```docker
docker run -it -d -v /path/to/local:/app bioinformatics_demo
docker run -it -d -v /Users/rjing/Desktop/test:/app bioinformatics_demo
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

```

## Memory Allocation

The recommended memory allocation for running this Docker image is at least 4GB. Depending on the size of the input data and the complexity of the analysis, you may need to allocate more memory.

To allocate more memory, use the `-m` flag when running the Docker container:

```docker
docker run --rm -m 8g xyz-analysis
```


This will allocate 8GB of memory to the container.

## Execution Time

The execution time of the analysis pipeline depends on various factors, including the size of the input data and the complexity of the analysis. It is difficult to provide an exact estimate without specific details.

However, as a rough estimate, the analysis pipeline typically takes around X hours to complete. Please note that this is just an estimate and the actual execution time may vary.

For any further questions or issues, please contact the XYZ Analysis team at example@example.com.
