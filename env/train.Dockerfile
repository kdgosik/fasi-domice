FROM registry.gitlab.com/fl84inc/comp-bio:main

RUN pip install --no-cache-dir \
    mord \
    torch \
    torchvision \
    torchaudio
