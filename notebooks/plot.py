import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import cmocean

import utils


long_name = dict(
    Cant='C$_{ant}$',
    pCFC11='pCFC-11',
    pCFC12='pCFC-12',    
    Del14C='$\Delta^{14}$C',
    ALK='Alkalinity',
    photoC_TOT_zint_100m='NPP',
    POC_FLUX_100m='POC export',
    POP_FLUX_100m='POP export',
    SiO2_FLUX_100m='Opal export',
    CaCO3_FLUX_100m='CaCO3 export',
)
units = dict(
    Cant='mmol m$^{-3}$',
    pCFC11='patm',
    pCFC12='patm',
    Del14C='‰',    
)


def za_obs_comparison(ds_zonal_mean, field, levels, levels_bias, col_name):
    """produce a plot of zonal mean fields from model, obs, and bias"""
    cmap_field = cmocean.cm.dense
    cmap_bias = cmocean.cm.balance
    
    extend = 'max' if levels[0] == 0. else 'both'
    contour_spec = {
        field: dict(
            levels=levels,
            extend=extend,
            cmap=cmap_field,
            norm=colors.BoundaryNorm(levels, ncolors=cmap_field.N),
        ),
        f'{field}_bias': dict(
            levels=levels_bias,
            extend='both',        
            cmap=cmap_bias,
            norm=colors.BoundaryNorm(levels_bias, ncolors=cmap_bias.N),
        ),
    }
    contour_spec.update({
        f'{v}_obs': contour_spec[v] 
        for v in contour_spec.keys() if '_bias' not in v
    })

    gss = []
    plot_basins = ['Atlantic Ocean', 'Pacific Ocean']
        
    fig = plt.figure(figsize=(11, 6))
    gs = gridspec.GridSpec(len(plot_basins), 3)

    axs_surf = []
    axs_deep = []
    cfs_surf = []
    for row, basin in enumerate(plot_basins):

        axs_surf_row = []
        axs_deep_row = []
        cfs_surf_row = []
        for col, v in enumerate([field, f'{field}_obs', f'{field}_bias']):
            gsi = gridspec.GridSpecFromSubplotSpec(100, 1, subplot_spec=gs[row, col])

            ax_surf = fig.add_subplot(gsi[:45, 0])
            ax_deep = fig.add_subplot(gsi[47:, 0])

            axs_surf_row.append(ax_surf)
            axs_deep_row.append(ax_deep)

            cf = ax_surf.contourf(
                ds_zonal_mean.lat_t, ds_zonal_mean.z_t * 1e-2, 
                ds_zonal_mean[v].sel(basins=basin),
                **contour_spec[v]
            )
            cfs_surf_row.append(cf)

            ax_deep.contourf(
                ds_zonal_mean.lat_t, ds_zonal_mean.z_t * 1e-2, 
                ds_zonal_mean[v].sel(basins=basin),
                **contour_spec[v]
            )

            ax_surf.contour(
                ds_zonal_mean.lat_t, ds_zonal_mean.z_t * 1e-2, 
                ds_zonal_mean[v].sel(basins=basin),
                levels=contour_spec[v]['levels'], colors='k', linewidths=0.2,
            )

            ax_deep.contour(
                ds_zonal_mean.lat_t, ds_zonal_mean.z_t * 1e-2, 
                ds_zonal_mean[v].sel(basins=basin),
                levels=contour_spec[v]['levels'], colors='k', linewidths=0.2,
            )


            ax_surf.set_ylim([1000., 0.])
            ax_surf.set_yticks(np.arange(0, 1000, 200))
            ax_surf.set_xticklabels([])
            ax_surf.xaxis.set_ticks_position('top')
            ax_surf.set_xticks(np.arange(-90, 110, 30))

            ax_deep.set_ylim([5000., 1000.])
            ax_deep.xaxis.set_ticks_position('bottom')
            ax_deep.set_xticks(np.arange(-90, 110, 30))

            if col == 0:
                ax_deep.set_ylabel('Depth [m]')
                ax_deep.yaxis.set_label_coords(-0.18, 1.05)
            else:
                ax_surf.set_yticklabels('')
                ax_deep.set_yticklabels('')

            if row == 1:
                ax_deep.set_xlabel('Latitude [$^\circ$N]')
            else:
                ax_surf.set_title(f'{long_name[field]} {col_name[col]}', loc='left')
                ax_surf.set_title(units[field], loc='right')                            
                ax_deep.set_xticklabels('')


        axs_surf.append(axs_surf_row)
        axs_deep.append(axs_deep_row)
        cfs_surf.append(cfs_surf_row)

    gs.update(left=0.11, right=0.89, wspace=0.08,hspace=0.075)

    #-- shift the right two columns over to make room for colorbar
    offset = 0.05
    for i in range(2):
        for j in range(2, 3):
            p0 = axs_surf[i][j].get_position()
            axs_surf[i][j].set_position([p0.x0+offset,p0.y0,p0.width,p0.height])

            p0 = axs_deep[i][j].get_position()
            axs_deep[i][j].set_position([p0.x0+offset,p0.y0,p0.width,p0.height])


    #-- add colorbars
    for i in range(2):
        for j in [1, 2]:
            p0 = axs_surf[i][j].get_position()
            p1 = axs_deep[i][j].get_position()

            cbaxes = fig.add_axes([p1.x0 + p1.width + 0.01, 
                                   p1.y0 + 0.0, 
                                   0.01, 
                                   p0.height + p1.height - 0.0])
            cb = fig.colorbar(cfs_surf[i][j], cax=cbaxes)

    fig.text(0.03, 0.6, plot_basins[0],
             fontsize=14.,
             fontweight = 'semibold',rotation=90);

    fig.text(0.03, 0.2, plot_basins[1],
             fontsize=14.,
             fontweight = 'semibold',rotation=90);

    utils.label_plots(fig, [ax for ax_row in axs_deep for ax in ax_row] , xoff=0.005, yoff=-0.18)


    
def nice_levels(da, max_steps=30, outside=False):
    """
    Return nice contour levels
    outside indicates whether the contour should be inside or outside the bounds
    """
    import sys
    
    cmin = da.min().values
    cmax = da.max().min()    
    
    table = [1., 2., 2.5, 4., 5., 10., 20., 25., 40., 50.,
             100., 200., 250., 400., 500.]
    npts = len(table)

    d = 10.**(np.floor(np.log10(cmax-cmin))-2.)

    u = sys.float_info.max
    step_size = sys.float_info.max
    am2 = 0.
    ax2 = 0.
    if outside:
        for i in range(npts):
            t = table[i] * d
            am1 = np.floor(cmin/t) * t
            ax1 = np.ceil(cmax/t) * t

            if (i == npts-1 and step_size == u) or ( (t <= step_size) and ((ax1-am1)/t <= (max_steps-1)) ):
                step_size = t
                ax2 = ax1
                am2 = am1
    else:
        for i in range(npts):
            t = table[i] * d
            am1 = np.ceil(cmin/t) * t
            ax1 = np.floor(cmax/t) * t

            if (i == npts-1 and step_size == u) or ( (t <= step_size) and ((ax1-am1)/t <= (max_steps-1)) ):
                step_size = t
                ax2 = ax1
                am2 = am1

    min_out = am2
    max_out = ax2

    return np.arange(min_out, max_out + step_size, step_size)