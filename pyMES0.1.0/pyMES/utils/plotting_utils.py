import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

class PlottingUtils:
    """Comprehensive plotting utility class for scientific visualization."""

    @staticmethod
    def plot_series(x, y, xlabel="X", ylabel="Y", title="Plot", grid=True, save_path=None, show=True, yerr=None):
        plt.figure()
        plt.errorbar(x, y, yerr=yerr, fmt='-o')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        if grid:
            plt.grid(True, which='both', ls='--', alpha=0.6)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()
        plt.close()

    @staticmethod
    def plot_multiple_series(x, ys, labels, xlabel="X", ylabel="Y", title="Multiple Series Plot", grid=True, save_path=None, show=True, yerrs=None):
        plt.figure()
        for idx, (y, label) in enumerate(zip(ys, labels)):
            err = yerrs[idx] if yerrs is not None else None
            plt.errorbar(x, y, yerr=err, marker='o', label=label)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        if grid:
            plt.grid(True, which='both', ls='--', alpha=0.6)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()
        plt.close()

    @staticmethod
    def plot_scatter(x, y, xlabel="X", ylabel="Y", title="Scatter Plot", grid=True, save_path=None, show=True, yerr=None, xerr=None):
        plt.figure()
        plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', linestyle='None')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        if grid:
            plt.grid(True, which='both', ls='--', alpha=0.6)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()
        plt.close()

    @staticmethod
    def plot_bar(x, heights, xlabel="X", ylabel="Y", title="Bar Plot", grid=True, save_path=None, show=True, yerr=None):
        plt.figure()
        plt.bar(x, heights, yerr=yerr, capsize=5)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        if grid:
            plt.grid(axis='y', linestyle='--', alpha=0.6)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()
        plt.close()

    @staticmethod
    def plot_double_y_bar(x, heights1, heights2, ylabel1, ylabel2, xlabel="X", title="Double Y-Axis Bar", grid=True, save_path=None, show=True, yerr1=None, yerr2=None):
        fig, ax1 = plt.subplots()
        width = 0.35
        idx = np.arange(len(x))
        ax1.bar(idx - width/2, heights1, width, yerr=yerr1, label=ylabel1, capsize=5)
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel1)
        ax2 = ax1.twinx()
        ax2.bar(idx + width/2, heights2, width, yerr=yerr2, color='orange', label=ylabel2, capsize=5)
        ax2.set_ylabel(ylabel2)
        plt.title(title)
        if grid:
            ax1.grid(True, axis='y', linestyle='--', alpha=0.6)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()
        plt.close()

    @staticmethod
    def plot_double_y_series(x, y1, y2, ylabel1, ylabel2, xlabel="X", title="Double Y-Axis Series", grid=True, save_path=None, show=True, yerr1=None, yerr2=None):
        fig, ax1 = plt.subplots()
        ax1.errorbar(x, y1, yerr=yerr1, marker='o', color='blue', label=ylabel1)
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel1, color='blue')
        ax2 = ax1.twinx()
        ax2.errorbar(x, y2, yerr=yerr2, marker='s', color='orange', label=ylabel2)
        ax2.set_ylabel(ylabel2, color='orange')
        plt.title(title)
        if grid:
            ax1.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()
        plt.close()

    @staticmethod
    def plot_contour(X, Y, Z, xlabel="X", ylabel="Y", title="Contour Plot", colorbar_label="Z", save_path=None, show=True):
        plt.figure()
        cp = plt.contourf(X, Y, Z, cmap='viridis')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        cbar = plt.colorbar(cp)
        cbar.set_label(colorbar_label)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()
        plt.close()

    @staticmethod
    def plot_IV_curve(voltage, current, xlabel="Voltage (V)", ylabel="Current (A)", title="I-V Curve", grid=True, save_path=None, show=True, yerr=None):
        PlottingUtils.plot_series(voltage, current, xlabel, ylabel, title, grid, save_path, show, yerr)

    @staticmethod
    def plot_tafel_plot(current_density, overpotential, xlabel="log(Current Density) [A/m²]", ylabel="Overpotential [V]", title="Tafel Plot", grid=True, save_path=None, show=True, yerr=None):
        log_current_density = np.log10(current_density)
        PlottingUtils.plot_series(log_current_density, overpotential, xlabel, ylabel, title, grid, save_path, show, yerr)

    @staticmethod
    def plot_3d_surface(X, Y, Z, xlabel="X", ylabel="Y", zlabel="Z", title="3D Surface", cmap='viridis', save_path=None, show=True):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(X, Y, Z, cmap=cmap, edgecolor='none')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        ax.set_title(title)
        fig.colorbar(surf, shrink=0.5, aspect=8, pad=0.1)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()
        plt.close()

    @staticmethod
    def plot_3d_scatter(x, y, z, xlabel="X", ylabel="Y", zlabel="Z", title="3D Scatter", color='b', save_path=None, show=True):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c=color, marker='o')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        ax.set_title(title)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()
        plt.close()

    @staticmethod
    def plot_3d_wireframe(X, Y, Z, xlabel="X", ylabel="Y", zlabel="Z", title="3D Wireframe", color='gray', save_path=None, show=True):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(X, Y, Z, color=color)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        ax.set_title(title)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()
        plt.close()

