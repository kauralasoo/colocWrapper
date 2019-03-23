FROM nfcore/base
LABEL authors="Kaur Alasoo" \
      description="Docker image containing all requirements for colocWrapper R package and pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/colocWrapper-1.0dev/bin:$PATH
