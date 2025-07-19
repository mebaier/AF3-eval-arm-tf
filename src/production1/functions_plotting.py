import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def create_scatter_plot_colour(df: pd.DataFrame, x_metric: str, y_metric: str, color_metric: str, 
                        title: str = '', cmap: str = 'viridis', alpha: float = 0.7, size: int = 50, ax=None) -> None:
    """Create a scatter plot with two metrics and color by a third metric.
    
    This function creates a scatter plot between two specified metrics with points colored
    by a third metric.
    
    Args:
        df (pd.DataFrame): DataFrame containing the metrics to plot
        x_metric (str): Name of the column to plot on x-axis
        y_metric (str): Name of the column to plot on y-axis
        color_metric (str): Name of the column to use for point colors
        title (str, optional): Custom title for the plot. Defaults to ''.
        cmap (str, optional): Matplotlib colormap name. Defaults to 'viridis'.
        alpha (float, optional): Transparency of points. Defaults to 0.7.
        size (int, optional): Point size. Defaults to 50.
        ax (matplotlib.axes.Axes, optional): Axes object to plot on. If None, creates new figure.
        
    Returns:
        None
    """
    if ax is None:
        plt.figure(figsize=(10, 8))
        ax = plt.gca()
        show_plot = True
    else:
        show_plot = False
    
    scatter = ax.scatter(
        x=df[x_metric], 
        y=df[y_metric], 
        c=df[color_metric],
        cmap=cmap,
        alpha=alpha,
        s=size,
        edgecolors='w'  # White edge to make points stand out
    )
    
    if title == '':
        title = f'{y_metric} vs {x_metric}\nColored by {color_metric}'

    # Add a color bar to show the scale of the color metric
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(color_metric, fontsize=10)

    # Add labels and title
    ax.set_xlabel(x_metric, fontsize=12)
    ax.set_ylabel(y_metric, fontsize=12)
    ax.set_title(title, fontsize=12)

    # Add a grid for better readability
    ax.grid(True, linestyle='--', alpha=0.7)

    # Add number of datapoints as a text box in the upper left corner
    num_points = len(df)
    ax.annotate(f'Datapoints: {num_points}', xy=(0.05, 0.95), xycoords='axes fraction', 
                fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

    # Show plot only if not using subplots
    if show_plot:
        plt.tight_layout()
        plt.show()
        
def create_scatter_plot(df: pd.DataFrame, x_metric: str, y_metric: str,
                        title: str = '', alpha: float = 0.7, size: int = 50, ax=None) -> None:
    """Create a scatter plot with two metrics.
    
    This function creates a scatter plot between two specified metrics.

    Args:
        df (pd.DataFrame): DataFrame containing the metrics to plot
        x_metric (str): Name of the column to plot on x-axis
        y_metric (str): Name of the column to plot on y-axis
        title (str, optional): Custom title for the plot. Defaults to ''.
        alpha (float, optional): Transparency of points. Defaults to 0.7.
        size (int, optional): Point size. Defaults to 50.
        ax (matplotlib.axes.Axes, optional): Axes object to plot on. If None, creates new figure.
        
    Returns:
        None
    """
    if ax is None:
        plt.figure(figsize=(10, 8))
        ax = plt.gca()
        show_plot = True
    else:
        show_plot = False
    
    scatter = ax.scatter(
        x=df[x_metric], 
        y=df[y_metric], 
        alpha=alpha,
        s=size,
        edgecolors='w'  # White edge to make points stand out
    )
    
    if title == '':
        title = f'{y_metric} vs {x_metric}'

    # Add labels and title
    ax.set_xlabel(x_metric, fontsize=12)
    ax.set_ylabel(y_metric, fontsize=12)
    ax.set_title(title, fontsize=12)

    # Add a grid for better readability
    ax.grid(True, linestyle='--', alpha=0.7)

    # Add number of datapoints as a text box in the upper left corner
    num_points = len(df)
    ax.annotate(f'Datapoints: {num_points}', xy=(0.05, 0.95), xycoords='axes fraction', 
                fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

    # Show plot only if not using subplots
    if show_plot:
        plt.tight_layout()
        plt.show()
