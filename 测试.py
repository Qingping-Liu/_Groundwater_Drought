# # encoding:utf-8
#
# import os
#
#
#
# # path = r"C:\Users\11072\Desktop\GRACE_dataset\huebei\111.gdb"
# # name = "NCP_P"
# # bound = path + "\\" + name
# #
# #
# # outpath = r"D:\GRACE_ini_datasets\initial_INPUT_data\_NCP\Provinces"
# #
# # # with arcpy.da.SearchCursor(bound, ["SHAPE@", "Pro"]) as rows:
# # #     for row in rows:
# # #         name = row[1] + ".shp"
# # #         print(name)
# # #         arcpy.FeatureClassToFeatureClass_conversion(row[0], outpath, name)
# # #         for
# #
# #
# # name= "hh"
# # print("%s" % name)
#
#
# import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
# from matplotlib import colors
#
# vegetables = ["cucumber", "tomato", "lettuce", "asparagus",
#               "potato", "wheat", "barley"]
# farmers = ["Farmer Joe", "Upland Bros.", "Smith Gardening",
#            "Agrifun", "Organiculture", "BioGoods Ltd.", "Cornylee Corp."]
#
# harvest = np.array([[0.8, 2.4, 2.5, 3.9, 0.0, 4.0, 0.0],
#                     [2.4, 0.0, 4.0, 1.0, 2.7, 0.0, 0.0],
#                     [1.1, 2.4, 0.8, 4.3, 1.9, 4.4, 0.0],
#                     [0.6, 0.0, 0.3, 0.0, 3.1, 0.0, 0.0],
#                     [0.7, 1.7, 0.6, 2.6, 2.2, 6.2, 0.0],
#                     [1.3, 1.2, 0.0, 0.0, 0.0, 3.2, 5.1],
#                     [0.1, 2.0, 0.0, 1.4, 0.0, 1.9, 6.3]])
#
# def heatmap(data, row_labels, col_labels, ax=None,
#             cbar_kw={}, cbarlabel="", **kwargs):
#     """
#     Create a heatmap from a numpy array and two lists of labels.
#
#     Parameters
#     ----------
#     data
#         A 2D numpy array of shape (M, N).
#     row_labels
#         A list or array of length M with the labels for the rows.
#     col_labels
#         A list or array of length N with the labels for the columns.
#     ax
#         A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
#         not provided, use current axes or create a new one.  Optional.
#     cbar_kw
#         A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
#     cbarlabel
#         The label for the colorbar.  Optional.
#     **kwargs
#         All other arguments are forwarded to `imshow`.
#     """
#     from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#
#     red = ListedColormap(["darkorange"])
#     blue= ListedColormap(["blue"])
#
#
#     if not ax:
#         ax = plt.gca()
#
#     bounds = np.array([0, 2.0, 2.5, 3.5, 10])
#     norm = colors.BoundaryNorm(boundaries=bounds, ncolors=4)
#     # Plot the heatmap
#     im = ax.imshow(data, **kwargs, norm=norm)
#
#     # Create colorbar
#     cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
#     cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
#
#     # Show all ticks and label them with the respective list entries.
#     ax.set_xticks(np.arange(data.shape[1]))
#     ax.set_yticks(np.arange(data.shape[0]))
#
#     # Let the horizontal axes labeling appear on top.
#     ax.tick_params(top=True, bottom=False,
#                    labeltop=True, labelbottom=False)
#
#     # Rotate the tick labels and set their alignment.
#     plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
#              rotation_mode="anchor")
#
#     # Turn spines off and create white grid.
#     ax.spines["right"].set_visible(False)
#     ax.spines["left"].set_visible(False)
#     ax.spines["top"].set_visible(False)
#     ax.spines["bottom"].set_visible(False)
#
#     ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
#     ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
#     ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
#     ax.tick_params(which="minor", bottom=False, left=False)
#
#     cbar = 0
#     return im, cbar
#
#
# def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
#                      textcolors=("black", "white"),
#                      threshold=None, **textkw):
#     """
#     A function to annotate a heatmap.
#
#     Parameters
#     ----------
#     im
#         The AxesImage to be labeled.
#     data
#         Data used to annotate.  If None, the image's data is used.  Optional.
#     valfmt
#         The format of the annotations inside the heatmap.  This should either
#         use the string format method, e.g. "$ {x:.2f}", or be a
#         `matplotlib.ticker.Formatter`.  Optional.
#     textcolors
#         A pair of colors.  The first is used for values below a threshold,
#         the second for those above.  Optional.
#     threshold
#         Value in data units according to which the colors from textcolors are
#         applied.  If None (the default) uses the middle of the colormap as
#         separation.  Optional.
#     **kwargs
#         All other arguments are forwarded to each call to `text` used to create
#         the text labels.
#     """
#
#     if not isinstance(data, (list, np.ndarray)):
#         data = im.get_array()
#
#     # Normalize the threshold to the images color range.
#     if threshold is not None:
#         threshold = im.norm(threshold)
#     else:
#         threshold = im.norm(data.max())/2.
#
#     # Set default alignment to center, but allow it to be
#     # overwritten by textkw.
#     kw = dict(horizontalalignment="center",
#               verticalalignment="center")
#     kw.update(textkw)
#
#     # Get the formatter in case a string is supplied
#     if isinstance(valfmt, str):
#         valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)
#
#     # Loop over the data and create a `Text` for each "pixel".
#     # Change the text's color depending on the data.
#     texts = []
#     for i in range(data.shape[0]):
#         for j in range(data.shape[1]):
#             kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
#             text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
#             texts.append(text)
#
#     return texts
#
#
# fig = plt.figure()
# ax = fig.add_subplot()
# im, cbar = heatmap(harvest, vegetables, farmers, ax=ax,
#                    cmap="Dark2", cbarlabel="harvest [t/year]")
# #texts = annotate_heatmap(im, valfmt="{x:.1f} t")
#
# fig.tight_layout()
# plt.show()
# os.system("pause")
#
#
#
#


# import numpy as np
# import random
#
# arr = np.arange(10).tolist()
# print(random.shuffle(arr))

import numpy as np
arr = np.arange(10)
np.random.shuffle(arr)
print(arr)
# 输出：[2 4 0 1 8 7 3 6 5 9]
