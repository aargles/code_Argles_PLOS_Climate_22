# -*- coding: utf-8 -*-
"""
Created on 22/04/2022

Defined here are secondary definitions that are used locally in “produce_mod.py”, these definitions are primarily to keep code tidy and save space. Custom matplotlib styles (e.g. custom.mplstyle) used are defined in the “style” subfolder.

"""

# Credits
__author__ = "Arthur Argles, Jonathan Moore, and Peter Cox"
__credits__ = ["Arthur Argles","Jonathan Moore","Peter Cox"]
__license__ = "CCBY4.0"
__maintainer__ = "Arthur Argles"
__email__ = "A.Argles2@exeter.ac.uk"


import os
import numpy as np
import matplotlib.pyplot as plt

def plot_style(style_choice='custom.mplstyle',global_style=False):
    """
    Function selects the matplotlib .mplstyle format either from the local directory …/routine/styles, or from the default list (https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html) when global_style is set to True.

    Parameters
    ----------
    style_choice : str, optional
        The chosen mplstyle filename. The default is 'custom.mplstyle'.
    global_style : str, optional
        If true select a template style as defined in the matplotlibrc folder. Otherwise select the local copy at …/routine/styles. The default is False.
        
    Returns
    -------
    None.

    """
    
    if global_style == True:
        
        plt.style.use(style_choice)      
         
    else:
        
        style_folder = os.path.join(os.getcwd(),os.path.join('routines',\
                                                             'styles'))
        plt.style.use(style_folder+'/'+style_choice)
        
    return


def logmidpoint(x1,x2):
    """
    Determine where the middle logarithmic value between two points.

    Parameters
    ----------
    x1 : float
        A point in the x dimension.
    x2 : float
        A point in the x dimension.

    Returns
    -------
    x_logmid : float
        The logarithmic midpoint between x1 and x2.

    """
    
    log_mid = 0.5*np.log(x1) + 0.5*np.log(x2)
    x_logmid = np.exp(log_mid)
    
    return x_logmid

from matplotlib import patches


def interpret_text_pos(text_pos_str,x,y,scale="linear"):
    """
    Function helps interpret the "text pos" input for figures from there data file.

    Parameters
    ----------
    text_pos_str : str
        Where to place a string relative to object. First word is vertical placement ("top", "center", or "bottom"), second word is horizontal placement ("left", "center", or "right), e.g. "top left" (there must be a space).
    x : array_like of size dimensions (2,)
        The minimum and maximum along the x-axis.
    y : array_like of dimensions = (2,)
        The minimum and maximum along the y-axis.
    scale : str, optional
        The scale of the axis: "linear" or "log". The default is "linear".

    Returns
    -------
    text_pos : array_like of dimensions = (2,)
        x and y coordinates of text annotation.
    text_ha : str
        horizontal alignment of annotation relative to text_pos
    text_va : TYPE
        vertical alignment of annotation relative to text_pos

    """
    
    if scale == 'linear':
        
        midpoint = (0.5*(x[0]+x[1]),0.5*(y[0]+y[1]))
    
    elif scale == 'log':
        
        midpoint = (logmidpoint(x[0],x[1]),\
                    logmidpoint(y[0],y[1]))        
        
    text_ha = 'center'
    text_va = 'bottom'            
    text_pos = [0,0]
    text_v_pos, text_h_pos = text_pos_str.split(' ')
    
    if text_v_pos == 'top':
        
        text_pos[1] = y[1]
                        
    elif text_v_pos == 'center':
        
        text_pos[1] = midpoint[1]
        text_va  = 'center'
        
    elif text_v_pos == 'bottom':
        
        text_pos[1] = y[0]
        text_va = 'top'
             
    if text_h_pos == 'left':
        
        text_pos[0] = x[0]
        text_ha = 'left'
        
    elif text_h_pos == 'center':
        
        text_pos[0] = midpoint[0]
        
    elif text_h_pos == 'right':
        
        text_pos[0] = x[1]
        text_ha = 'right'
        
    return text_pos, text_ha, text_va