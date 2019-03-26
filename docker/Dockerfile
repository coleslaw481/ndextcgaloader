# Use an official Python runtime as a parent image
FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get install -y build-essential 

RUN mkdir /tmp/ndextcgaloader
COPY ./ /tmp/ndextcgaloader/
RUN pip install /tmp/ndextcgaloader

RUN rm -rf /tmp/ndextcgaloader

RUN mkdir /ndextcgaloader
WORKDIR /ndextcgaloader

CMD /bin/echo "To use run: ndexloadtcga.py"
