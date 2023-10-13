FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/bioinf-tools/:/bioinf-tools/cortex/bin/:/bioinf-tools/cortex/scripts/analyse_variants/seq-align/bin/:/bioinf-tools/vcftools-0.1.15/install/bin:/bioinf-tools/enaBrowserTools/python3:/clockwork/scripts/:$PATH
ENV PERL5LIB=/bioinf-tools/vcftools-0.1.15/install/share/perl/5.30.0/:/bioinf-tools/cortex/scripts/analyse_variants/bioinf-perl/lib:/bioinf-tools/cortex/scripts/calling/:$PERL5LIB
ENV LANG=C.UTF-8
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
ENV NXF_VER=22.10.0

COPY scripts/install_dependencies.sh /install_dependencies.sh
RUN /install_dependencies.sh /bioinf-tools

ARG CLOCKWORK_DIR=/clockwork
COPY . $CLOCKWORK_DIR
RUN service mysql start \
  && service mysql status \
  && cd $CLOCKWORK_DIR/python \
  && python3 setup.py test \
  && python3 setup.py install \
  && service mysql stop \
  && apt remove -y mysql-server

CMD clockwork
