FROM granatumx/gbox-py-sdk:1.0.0

RUN pip install seaborn

COPY . .

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
