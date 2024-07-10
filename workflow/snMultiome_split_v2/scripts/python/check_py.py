import matplotlib
matplotlib.use('agg')
import numpy as np
import scrublet as scr
import scipy.io
import os
from scipy.sparse import csc_matrix
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from tabulate import tabulate
import h5py
import pandas as pd
from scipy.sparse import csc_matrix
import math
from sklearn.mixture import GaussianMixture
from scipy.optimize import root_scalar
from scipy.stats import norm
from collections import Counter
import random
import argparse
import subprocess
import sys
import glob
import snakemake

