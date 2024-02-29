#!/usr/bin/python3
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import numpy as np

def grouped_stacked_bar2(lists1, lists2, x_tick_labels, labels1, labels2, title, xlabel, ylabel, colors1, colors2, min_y, max_y, 
                optional_annot, path, custom_legend=False, **kwargs):
    """
    Create a stacked bar chart from 2 lists of lists ('lists1 and 2').
    """
    rc = {'ytick.labelsize': 30, 'xtick.labelsize': 35}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(25, 8))
    df1 = pd.DataFrame(lists1, index=x_tick_labels, columns=labels1)
    df2 = pd.DataFrame(lists2, index=x_tick_labels, columns=labels2)
    print(df1, df2)
    df1.plot.bar(stacked=True, ax=ax, color=colors1, position=1, **kwargs)
    df2.plot.bar(stacked=True, ax=ax, color=colors2, position=0, **kwargs)
    ax.set_xticklabels(x_tick_labels, rotation=90)
    
    # Add optional annotation above bars (ex: number of sno in each bar)
    n = len(labels1)
    
    def chunker(list_, size):
        # Get chunks of n elements from list_ where n=size
        return (list_[pos:pos + size] for pos in range(0, len(list_), size))
    
    # Get the cumulative height of each bar until the before last stacked bar
    previous_heigths = [0] * len(x_tick_labels)
    for i, chunks_index in enumerate(chunker(list(range(0, n*len(x_tick_labels))[:-len(x_tick_labels)]), len(x_tick_labels))):  # from 0 to the before last bars
        for index_ in chunks_index:
            bar = list(ax.patches)[index_]
            previous_heigths[index_ - len(x_tick_labels) * i] += bar.get_height()
    
    # Add the cumulative height of previous bars to the last stacked bar of each stack
    last_bars = [bar_ for j, bar_ in enumerate(ax.patches) if j in list(range(0, n*len(x_tick_labels))[-len(x_tick_labels):])]  # last (i.e. topmost) bars
    longest_bar = 0
    for i, bar in enumerate(last_bars):
        if bar.get_height() + previous_heigths[i] + 0.35 >= longest_bar:
            longest_bar = bar.get_height() + previous_heigths[i] + 0.35
            ax.text((bar.get_x()+0.35), 
                    (bar.get_height() + previous_heigths[i] + 0.35), optional_annot[i], fontsize=25, weight='bold', horizontalalignment='center')
        else:
            ax.text((bar.get_x()+0.35), 
                    (longest_bar), optional_annot[i], fontsize=25, weight='bold', horizontalalignment='center')
        
    if custom_legend == True:
        legend_list = []
        n_crit = len(set(labels1))
        for i, crit in enumerate(labels2):
            legend_element = mpatches.Patch(color=colors2[i], label=crit)
            legend_list.append(legend_element)
        for i, crit in enumerate(labels1[0:n_crit]):
            legend_element = mpatches.Patch(color=colors1[i], label=crit)
            legend_list.append(legend_element)
        plt.legend(handles=legend_list, loc=5, bbox_to_anchor=(0.75, 1.5), reverse=True, fontsize=25)
    else:
        plt.legend(fontsize=25, loc=5, bbox_to_anchor=(0.75, 1.5), reverse=True)
    plt.title(title, fontsize=40)
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=40)
    plt.autoscale()
    plt.margins(0.02)
    plt.ylim(min_y, max_y)
    plt.savefig(path, bbox_inches='tight', dpi=600)


def pie_multiple(y, x, count_list, labels, colors, ax_title, title, legend_title, path, **kwargs):
    """
    Creates x*y pie charts from a list of list (count_list) where each global
    element corresponds to a list of local elements (ex: percent per rank across
    8 model intersection).
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, axes = plt.subplots(y, x, figsize=(25, 8))
    plt.subplots_adjust(hspace=0.5)
    ax = axes.flatten()
    for i, element in enumerate(count_list):
        count_per_element = count_list[i][:]
        if count_per_element != [0, 0, 0]:
            ax[i].set_title(ax_title[i], fontdict={'fontsize': 12}, x=0.5, y=0.8)
            ax[i].pie(count_per_element, colors=colors, textprops={'fontsize': 19},
                        **kwargs)
            ax[i].axis('equal')
            white_circle = plt.Circle((0, 0), 0.4, color='white') #to create a donut chart
            ax[i].add_artist(white_circle) #to create a donut chart

    fig.suptitle(title, x=0.5, y=0.9, fontsize=25)
    fig.legend(labels=labels, loc='upper right', bbox_to_anchor=(0.9, 0.4),
                prop={'size': 20}, title=legend_title)
    plt.savefig(path, dpi=600)


def percent_count(count_list):
    """
    Create from a list of lists a percentage list of lists.
    """
    percent_list = []

    for i, cat in enumerate(count_list):
        temp = []
        for j, crit in enumerate(cat):
            total = sum(cat)
            if total != 0:
                percent = crit/total * 100
                percent = round(percent, 2)
            else:
                percent = 0
            temp.append(percent)
        percent_list.append(temp)

    return percent_list
