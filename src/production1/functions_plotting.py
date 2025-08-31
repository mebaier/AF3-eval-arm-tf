import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn
from typing import Set, Tuple, List

def create_scatter_plot_colour(df: pd.DataFrame, x_metric: str, y_metric: str, color_metric: str, 
                        title: str = '', cmap: str = 'viridis', alpha: float = 0.7, size: int = 50, ax=None, corr=False) -> None:
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
        corr (bool, optional): Whether to calculate and display Pearson correlation coefficient. Defaults to False.
        
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
    annotation_text = f'Datapoints: {num_points}'
    
    # Calculate and add correlation coefficient if requested
    if corr:
        # Filter out rows with missing values for correlation calculation
        valid_corr_data = df[[x_metric, y_metric]].dropna()
        if len(valid_corr_data) > 1:  # Need at least 2 points for correlation
            correlation = valid_corr_data[x_metric].corr(valid_corr_data[y_metric])
            annotation_text += f'\nPearson r: {correlation:.3f}'
        else:
            annotation_text += f'\nPearson r: N/A (insufficient data)'
    
    ax.annotate(annotation_text, xy=(0.05, 0.95), xycoords='axes fraction', 
                fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

    # Show plot only if not using subplots
    if show_plot:
        plt.tight_layout()
        plt.show()
        
def create_scatter_plot(df: pd.DataFrame, x_metric: str, y_metric: str,
                        title: str = '', alpha: float = 0.7, size: int = 50, ax=None, corr=False) -> None:
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
        corr (bool, optional): Whether to calculate and display Pearson correlation coefficient. Defaults to False.
        
    Returns:
        None
    """
    if ax is None:
        plt.figure(figsize=(10, 8))
        ax = plt.gca()
        show_plot = True
    else:
        show_plot = False
    
    # Filter out rows with missing values in either metric
    valid_data = df[[x_metric, y_metric]].dropna()
    
    # Check if we have any valid data left
    if len(valid_data) == 0:
        ax.text(0.5, 0.5, f'No valid data points found.\nAll values missing for {x_metric} or {y_metric}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=12, 
                bbox=dict(boxstyle="round,pad=0.5", fc="yellow", alpha=0.7))
        if title == '':
            title = f'{y_metric} vs {x_metric} (No Data)'
        ax.set_title(title, fontsize=12)
        if show_plot:
            plt.tight_layout()
            plt.show()
        return
    
    scatter = ax.scatter(
        x=valid_data[x_metric], 
        y=valid_data[y_metric], 
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

    # Add number of valid datapoints as a text box in the upper left corner
    num_valid_points = len(valid_data)
    num_total_points = len(df)
    if num_valid_points == num_total_points:
        annotation_text = f'Datapoints: {num_valid_points}'
    else:
        annotation_text = f'Datapoints: {num_valid_points}/{num_total_points} (valid/total)'
    
    # Calculate and add correlation coefficient if requested
    if corr:
        correlation = valid_data[x_metric].corr(valid_data[y_metric])
        annotation_text += f'\nPearson r: {correlation:.3f}'
    
    ax.annotate(annotation_text, xy=(0.05, 0.95), xycoords='axes fraction', 
                fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

    # Show plot only if not using subplots
    if show_plot:
        plt.tight_layout()
        plt.show()

def create_heatmap(df: pd.DataFrame, x_metric: str, y_metric: str,
                   title: str = '', bins: int | list = 30, cmap: str = 'viridis', ax=None, corr=False, scatter_threshold: int = 2, bin_size: float | None = None) -> None:
    """Create a heatmap (2D histogram) with two metrics.
    
    This function creates a heatmap between two specified metrics showing the density
    of data points in different regions. Points in bins with fewer than scatter_threshold
    counts are overlaid as individual scatter points.

    Args:
        df (pd.DataFrame): DataFrame containing the metrics to plot
        x_metric (str): Name of the column to plot on x-axis
        y_metric (str): Name of the column to plot on y-axis
        title (str, optional): Custom title for the plot. Defaults to ''.
        bins (int, optional): Number of bins for the heatmap. Defaults to 30.
        cmap (str, optional): Matplotlib colormap name. Defaults to 'viridis'.
        ax (matplotlib.axes.Axes, optional): Axes object to plot on. If None, creates new figure.
        corr (bool, optional): Whether to calculate and display Pearson correlation coefficient. Defaults to False.
        scatter_threshold (int, optional): Minimum number of points in a bin to show as heatmap only. 
                                         Bins with fewer points show individual scatter points. Defaults to 10.
        bin_size (float, optional): Size of each bin. If provided, overrides the bins parameter by calculating 
                                   the number of bins based on data range. Defaults to None.
        
    Returns:
        None
    """
    if ax is None:
        plt.figure(figsize=(10, 8))
        ax = plt.gca()
        show_plot = True
    else:
        show_plot = False
    
    # Filter out rows with missing values in either metric
    valid_data = df[[x_metric, y_metric]].dropna()
    
    # Check if we have any valid data left
    if len(valid_data) == 0:
        ax.text(0.5, 0.5, f'No valid data points found.\nAll values missing for {x_metric} or {y_metric}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=12, 
                bbox=dict(boxstyle="round,pad=0.5", fc="yellow", alpha=0.7))
        if title == '':
            title = f'{y_metric} vs {x_metric} Heatmap (No Data)'
        ax.set_title(title, fontsize=12)
        if show_plot:
            plt.tight_layout()
            plt.show()
        return
    
    # Calculate bins based on bin_size if provided
    if bin_size is not None:
        x_range = valid_data[x_metric].max() - valid_data[x_metric].min()
        y_range = valid_data[y_metric].max() - valid_data[y_metric].min()
        bins_x = max(1, int(x_range / bin_size))
        bins_y = max(1, int(y_range / bin_size))
        bins = [bins_x, bins_y]
    
    # Create 2D histogram (heatmap)
    counts, xedges, yedges, im = ax.hist2d(
        valid_data[x_metric], 
        valid_data[y_metric], 
        bins=bins,
        cmap=cmap
    )
    
    # Find points in bins with less than 10 counts and overlay as scatter points
    # Determine which bin each point belongs to
    x_indices = np.digitize(valid_data[x_metric], xedges) - 1
    y_indices = np.digitize(valid_data[y_metric], yedges) - 1
    
    # Clip indices to valid range
    x_indices = np.clip(x_indices, 0, len(xedges) - 2)
    y_indices = np.clip(y_indices, 0, len(yedges) - 2)
    
    # Create mask for points in low-density bins (< scatter_threshold points)
    low_density_mask = counts[x_indices, y_indices] < scatter_threshold
    
    # Overlay scatter points for low-density regions
    if np.any(low_density_mask):
        low_density_data = valid_data[low_density_mask]
        ax.scatter(
            low_density_data[x_metric],
            low_density_data[y_metric],
            c='white',
            s=30,
            alpha=0.8,
            edgecolors='black',
            linewidth=0.5,
            zorder=5  # Ensure points are on top
        )
    
    if title == '':
        title = f'{y_metric} vs {x_metric} Heatmap + Scatter'

    # Add colorbar to show the count scale
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Count', fontsize=10)

    # Add labels and title
    ax.set_xlabel(x_metric, fontsize=12)
    ax.set_ylabel(y_metric, fontsize=12)
    ax.set_title(title, fontsize=12)

    # Add a grid for better readability
    ax.grid(True, linestyle='--', alpha=0.3)

    # Add number of valid datapoints as a text box in the upper left corner
    num_valid_points = len(valid_data)
    num_total_points = len(df)
    if num_valid_points == num_total_points:
        annotation_text = f'Datapoints: {num_valid_points}'
    else:
        annotation_text = f'Datapoints: {num_valid_points}/{num_total_points} (valid/total)'
    
    # Calculate and add correlation coefficient if requested
    if corr:
        correlation = valid_data[x_metric].corr(valid_data[y_metric])
        annotation_text += f'\nPearson r: {correlation:.3f}'
    
    ax.annotate(annotation_text, xy=(0.05, 0.95), xycoords='axes fraction', 
                fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

    # Show plot only if not using subplots
    if show_plot:
        plt.tight_layout()
        plt.show()

def create_overlaid_AF_violin_plot(dataset1, dataset1_name, dataset2, dataset2_name, 
                               figsize=(12, 6)):
    """
    Create overlaid violin plots comparing AlphaFold metrics between two datasets.
    
    Parameters:
    -----------
    dataset1 : pd.DataFrame
        First dataset containing the metrics
    dataset1_name : str
        Label for the first dataset
    dataset2 : pd.DataFrame  
        Second dataset containing the metrics
    dataset2_name : str
        Label for the second dataset
    figsize : tuple
        Figure size as (width, height)
        
    Returns:
    --------
    fig, ax : matplotlib objects
        Figure and axes objects for further customization
    """
    # Hard-coded metrics
    metrics = ['ranking_score', 'iptm', 'ptm']
    
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
    melted_data1 = pd.melt(dataset1[metrics], 
                           var_name='metric', value_name='value')
    melted_data1['dataset'] = dataset1_name
    
    melted_data2 = pd.melt(dataset2[metrics], 
                           var_name='metric', value_name='value')
    melted_data2['dataset'] = dataset2_name
    
    # Combine datasets
    combined_data = pd.concat([melted_data1, melted_data2], ignore_index=True)
    
    # Create overlaid violin plots
    sns.violinplot(data=combined_data, x='metric', y='value', hue='dataset', 
                   split=True, ax=ax, palette=['lightblue', 'lightcoral'])
    
    ax.set_title('Distribution of AlphaFold Metrics - Dataset Comparison')
    ax.set_ylabel('Metric Values')
    ax.set_xlabel('AlphaFold Metrics')
    
    # Add text box with number of data points for both datasets
    n_points1 = len(dataset1)
    n_points2 = len(dataset2)
    text_content = f'{dataset1_name}: n = {n_points1}\n{dataset2_name}: n = {n_points2}'
    ax.text(0.02, 0.98, text_content, transform=ax.transAxes, 
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            verticalalignment='top', fontsize=10)
    
    plt.tight_layout()
    
    return fig, ax

def create_single_AF_violin_plot(dataset, dataset_name, figsize=(10, 6)):
    """
    Create a single violin plot for AlphaFold metrics from one dataset.
    
    Parameters:
    -----------
    dataset : pd.DataFrame
        Dataset containing the metrics
    dataset_name : str
        Label for the dataset
    figsize : tuple
        Figure size as (width, height)
   
    Returns:
    --------
    fig, ax : matplotlib objects
        Figure and axes objects for further customization
    """
    # Hard-coded metrics
    metrics = ['ranking_score', 'iptm', 'ptm']
    
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
    melted_data = pd.melt(dataset[metrics], 
                          var_name='metric', value_name='value')
    
    # Create violin plot
    sns.violinplot(data=melted_data, x='metric', y='value', ax=ax)
    ax.set_title(f'Distribution of AlphaFold Metrics - {dataset_name}')
    ax.set_ylabel('Metric Values')
    ax.set_xlabel('AlphaFold Metrics')
    
    # Add text box with number of data points
    n_points = len(dataset)
    ax.text(0.02, 0.98, f'n = {n_points}', transform=ax.transAxes, 
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            verticalalignment='top', fontsize=10)
    
    plt.tight_layout()
    
    return fig, ax

def create_violin_plot(dataset, metric, title='', xlabel='', ylabel='', figsize=(10, 6)):
    """
    Create a single violin plot for a specific metric from one dataset.
    
    Parameters:
    -----------
    dataset : pd.DataFrame
        Dataset containing the metric
    metric : str
        Name of the metric column to plot
    title : str
        Title for the plot
    xlabel : str
        Label for x-axis
    ylabel : str
        Label for y-axis
    figsize : tuple
        Figure size as (width, height)
        
    Returns:
    --------
    fig, ax : matplotlib objects
        Figure and axes objects for further customization
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
    melted_data = pd.melt(dataset[[metric]], 
                          var_name='metric', value_name='value')
    
    # Create violin plot
    sns.violinplot(data=melted_data, x='metric', y='value', ax=ax)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    # Make x-axis tighter by setting limits
    ax.set_xlim(-0.5, 0.5)
    
    # Add text box with number of data points
    n_points = len(dataset)
    ax.text(0.02, 0.98, f'n = {n_points}', transform=ax.transAxes, 
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            verticalalignment='top', fontsize=10)
    
    plt.tight_layout()
    
    return fig, ax

def create_double_violin_plot(dataset1, dataset1_name, dataset2, dataset2_name, 
                              metric1, metric2, figsize=(12, 6)):
    """
    Create a split violin plot comparing two specific metrics between two datasets.
    Shows dataset1 with metric1 on the left side and dataset2 with metric2 on the right side of one violin.
    
    Parameters:
    -----------
    dataset1 : pd.DataFrame
        First dataset containing the metrics
    dataset1_name : str
        Label for the first dataset
    dataset2 : pd.DataFrame  
        Second dataset containing the metrics
    dataset2_name : str
        Label for the second dataset
    metric1 : str
        Name of the first metric column to plot (from dataset1)
    metric2 : str
        Name of the second metric column to plot (from dataset2)
    figsize : tuple
        Figure size as (width, height)
        
    Returns:
    --------
    fig, ax : matplotlib objects
        Figure and axes objects for further customization
    """
    fig = plt.figure(figsize=figsize)
    ax = plt.gca()
    
    # Create data for dataset1 with metric1
    data1 = pd.DataFrame({
        'metric': 'comparison',
        'value': dataset1[metric1],
        'dataset': f'{dataset1_name}_{metric1}'
    })
    
    # Create data for dataset2 with metric2  
    data2 = pd.DataFrame({
        'metric': 'comparison',
        'value': dataset2[metric2],
        'dataset': f'{dataset2_name}_{metric2}'
    })
    
    # Combine datasets
    combined_data = pd.concat([data1, data2], ignore_index=True)
    
    # Create split violin plot
    sns.violinplot(data=combined_data, x='metric', y='value', hue='dataset',
                   split=True, ax=ax, palette=['lightblue', 'lightcoral'])
    
    ax.set_title(f'Distribution Comparison: {dataset1_name} ({metric1}) vs {dataset2_name} ({metric2})')
    ax.set_ylabel('Metric Values')
    ax.set_xlabel('')
    ax.set_xticklabels([metric1]) # FIXME: what if metric1 != metric2
    
    # Add text box with number of data points for both datasets
    n_points1 = len(dataset1)
    n_points2 = len(dataset2)
    text_content = f'{dataset1_name}: n = {n_points1}\n{dataset2_name}: n = {n_points2}'
    ax.text(0.02, 0.98, text_content, transform=ax.transAxes, 
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            verticalalignment='top', fontsize=10)
    
    plt.tight_layout()
    
    return fig, ax

def plot_roc(df: pd.DataFrame, pos: Set[str], neg: Set[str], param_name: str, min_val: float, max_val: float, 
             sampling: int = 1000, direction: str = 'up', plot_title: str = '', ax=None) -> Tuple[List[float], List[float]]:
    """Calculate and plot the ROC curve based on a specified parameter.
    
    This function evaluates the performance of a binary classifier by varying a threshold parameter
    and calculating the True Positive Rate (TPR) and False Positive Rate (FPR) at each threshold.
    
    Args:
        df (pd.DataFrame): DataFrame containing the parameter to evaluate and pair_id column
        pos (Set[str]): Set of positive example pair_ids (ground truth positive cases)
        neg (Set[str]): Set of negative example pair_ids (ground truth negative cases)
        param_name (str): Name of the column in df to use as the classification parameter
        min_val (float): Minimum threshold value to evaluate
        max_val (float): Maximum threshold value to evaluate
        sampling (int, optional): Number of threshold points to sample between min and max. Defaults to 1000.
        direction (str, optional): Direction of classification - 'up' means values >= threshold are positive,
                                  'down' means values < threshold are positive. Defaults to 'up'.
        plot_title (str, optional): Title for the ROC curve plot. If empty, a default title is used. Defaults to ''.
        ax (matplotlib.axes.Axes, optional): Axes object to plot on. If None, creates new figure.
                                   
    Returns:
        Tuple[List[float], List[float]]: Lists of FPR and TPR values that make up the ROC curve
    """
    if ax is None:
        # Print the size of positive and negative sets for verification    
        print(f"Number of positive examples: {len(pos)}")
        print(f"Number of negative examples: {len(neg)}")
        plt.figure(figsize=(8, 8))
        ax = plt.gca()
        show_plot = True
    else:
        show_plot = False

    # Calculate step size based on range and sampling
    step = (max_val - min_val) / sampling
    if step <= 0:
        raise ValueError("max_val must be greater than min_val")
    
    # Lists to store TPR and FPR values
    tpr_list = []
    fpr_list = []
    
    for i in range(sampling + 1):
        sep_val = min_val + (i * step)
        
        if direction == 'up':
            calc_pos = set(df[df[param_name] >= sep_val]['Entry ID'].tolist())
            calc_neg = set(df[df[param_name] < sep_val]['Entry ID'].tolist())
        elif direction == 'down':
            calc_pos = set(df[df[param_name] < sep_val]['Entry ID'].tolist())
            calc_neg = set(df[df[param_name] >= sep_val]['Entry ID'].tolist())
        else:
            raise ValueError("direction must be either 'up' or 'down'")
        
        # Calculate confusion matrix values
        TP_num = len(calc_pos.intersection(pos))
        FP_num = len(calc_pos.intersection(neg))
        TN_num = len(calc_neg.intersection(neg))
        FN_num = len(calc_neg.intersection(pos))
        
        # Avoid division by zero
        TPR = TP_num / max(TP_num + FN_num, 1)
        FPR = FP_num / max(FP_num + TN_num, 1)
        
        tpr_list.append(TPR)
        fpr_list.append(FPR)
    
    # Plot the ROC curve
    ax.plot(fpr_list, tpr_list, 'b-', linewidth=2)
    ax.plot([0, 1], [0, 1], 'k--', linewidth=2)  # Diagonal line representing random guess
    
    auc_value = sklearn.metrics.auc(fpr_list, tpr_list)
    
    # Add labels and title
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    
    if plot_title:
        title = plot_title
    else:
        title = f'ROC Curve for {param_name}'
    
    ax.set_title(f'{title}\nAUC = {auc_value:.3f}', fontsize=12)
    
    # Add grid and improve appearance
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    
    # Show plot only if not using subplots
    if show_plot:
        plt.tight_layout()
        plt.show()
    
    return fpr_list, tpr_list