## CRISPRtile

## Git Clone
    git clone https://github.com/jasoncngo/CRISPRtile.git
    cd CRISPRtile

## Mac Installation
    conda create -n CRISPRtile-venv python=3.10 python-kaleido=0.1.* -y
    conda activate CRISPRtile-venv
    conda install -c conda-forge mamba -y 
    mamba install -c conda-forge autogluon.tabular=0.7.0 mdtraj=1.9.7 -y 

## Windows Installation
    conda create -n CRISPRtile-venv python=3.10 python-kaleido=0.1.* mdtraj=1.9.7 cudatoolkit=11.3 -y 
    conda activate CRISPRtile-venv
    pip install -U pip
    pip install -U setuptools wheel
    pip install autogluon.tabular[all]==0.7.0 torch==1.13.1+cu116 torchvision==0.14.1+cu116 --extra-index-url https://download.pytorch.org/whl/cu116

## Contact
If you have any questions, please contact me at: <a href="mailto:jcn2144@columbia.edu">jcn2144@columbia.edu</a>
