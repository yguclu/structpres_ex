# structpres_ex
Exercises for lectures on [Structure Preserving Methods on Staggered Grids](http://www-m16.ma.tum.de/Allgemeines/StructPresStag19) (TUM, summer semester 2019).


## Lecture format

For each exercise assignment, the solution will be provided in the form of an IPython notebook (a file with extension *.ipynb*) that mixes

- [Markdown](https://www.markdownguide.org/getting-started) formatted text
- [LaTeX](https://www.latex-project.org/) formulas
- [TikZ](https://github.com/pgf-tikz/pgf) diagrams
- [Python3](https://www.python.org/about/) code
- [IPython](https://ipython.org/) "magic" commands

The notebook is used for interactive data exploration, and it provides a high-level interface to dedicated Python3 modules and scripts (which take care of the real number crunching). In fact, it is also possible to run the scripts from a system terminal or an IPython console.

In order to effectively navigate through a solution notebook, it is *strongly* recommended to open it with [JupyterLab](https://jupyterlab.readthedocs.io/en/stable/user/interface.html). JupyterLab is included in the latest Python3 [Anaconda distribution](https://www.anaconda.com/distribution/), or it may be installed using [pip](https://pip.pypa.io/en/stable/).

The JupyterLab distribution ships with a lot of built-in functionality, but we will use a few additional extensions to improve our workflow. The next section explains how to set up your environment.


## JupyterLab environment setup

Install JupyterLab if you do not have it

    $ pip3 install --user jupyterlab

Install the `ipympl` module

    $ pip3 install --user ipympl

Install `nodejs` and `npm` libraries with one of the two following options

- Download binary from https://nodejs.org/en/download/, or
- Use your Linux package manager, e.g.
     
      $ sudo apt install nodejs
      $ sudo apt install npm

Install the Matplotlib widget extension

    $ jupyter labextension install @jupyter-widgets/jupyterlab-manager
    $ jupyter labextension install jupyter-matplotlib

**(Optional)** Install the variable inspector extension

    $ jupyter labextension install @lckr/jupyterlab_variableinspector

**(Optional)** Install the table of contents (TOC) extension

    $ jupyter labextension install @jupyterlab/toc

**(Optional)** Install the IPython magics for generating figures with TikZ

    $ pip3 install git+git://github.com/mkrphys/ipython-tikzmagic.git
