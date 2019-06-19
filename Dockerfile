# Use an official Python runtime as a parent image
FROM ubuntu

# Install any needed packages specified in requirements.txt
RUN apt-get update && apt-get install -y git cmake g++ libgmp3-dev

# Set the working directory to /app
WORKDIR /app

# Download and compile TriWild
RUN git clone https://github.com/wildmeshing/TriWild/
WORKDIR /app/TriWild/build
RUN cmake .. && make

WORKDIR /data

#ENTRYPOINT ["/app/TriWild/build/TriWild"]

## Create TriWild image with:
# docker build -t triwild .
## Run TriWild with:
# docker run --rm -v "$(pwd)":/data triwild [TriWild arguments]
