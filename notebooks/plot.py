from itertools import product

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import cmocean

import utils

# this is really lame
long_name = dict(
    Cant='C$_{ant}$',
    Cant_v1='C$_{ant}$',
    Cant_v1pGruber2019='C$_{ant}$',
    pCFC11='pCFC-11',
    pCFC12='pCFC-12',    
    Del14C='$\Delta^{14}$C',
    ALK='Alkalinity',
    photoC_TOT_zint_100m='NPP',
    POC_FLUX_100m='POC export',
    POP_FLUX_100m='POP export',
    SiO2_FLUX_100m='Opal export',
    CaCO3_FLUX_100m='CaCO3 export',
    sDIC='sDIC',
    sALK='sAlk',    
)
units = dict(
    Cant='mmol m$^{-3}$',
    Cant_v1='mmol m$^{-3}$',
    Cant_v1pGruber2019='mmol m$^{-3}$',    
    pCFC11='patm',
    pCFC12='patm',
    Del14C='‰',
    sDIC='mmol m$^{-3}$',
    sALK='meq m$^{-3}$',
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


    
def nice_levels(da, max_steps=30, outside=False, percentile=[0., 100.]):
    """
    Return nice contour levels
    outside indicates whether the contour should be inside or outside the bounds
    """
    import sys
    
    if percentile[0] == 0.:
        cmin = da.min()
    else:
        cmin = np.nanpercentile(da, percentile[0])
    
    if percentile[1] == 100.:
        cmax = da.max()
    else:
        cmax = np.nanpercentile(da, percentile[1])    

    
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


def canvas(*args, figsize=(6, 4), use_gridspec=False, **gridspec_kwargs):

    assert len(args), 'Args required'
    assert len(args) <= 2, 'Too many args'
    
    if len(args) == 2:
        nrow = args[0]
        ncol = args[1]
    else:
        npanel = args[0]
        nrow = int(np.sqrt(npanel))
        ncol = int(npanel/nrow) + min(1, npanel%nrow)
    
    if use_gridspec:
        fig = plt.figure(figsize=(figsize[0]*ncol, figsize[1]*nrow)) #dpi=300)
        gs = gridspec.GridSpec(nrows=nrow, ncols=ncol, **gridspec_kwargs)
        axs = np.empty((nrow, ncol)).astype(object)
        for i, j in product(range(nrow), range(ncol)):
            axs[i, j] = plt.subplot(gs[i, j])
        return fig, axs
    else:
        return plt.subplots(
            nrow, ncol, 
            figsize=(figsize[0]*ncol, figsize[1]*nrow),                       
            constrained_layout=False,
            squeeze=False,
        )  
    
def moc(MOC):
    fig = plt.figure(figsize=(14, 4))
    gs = gridspec.GridSpec(
        1, len(MOC.transport_reg), 
        wspace=0.1,
    )

    lo, hi, dc = -40., 40., 4
    cnlevels = np.arange(lo, hi+dc, dc)

    axs_surf = []
    axs_deep = []
    for n, transport_reg in enumerate(MOC.transport_reg.values):
        gsi = gridspec.GridSpecFromSubplotSpec(100, 1, subplot_spec=gs[0, n])

        ax_surf = fig.add_subplot(gsi[:45, 0])
        ax_deep = fig.add_subplot(gsi[47:, 0])    

        ax_surf.set_facecolor('gainsboro')
        ax_deep.set_facecolor('gainsboro')

        axs_surf.append(ax_surf)
        axs_deep.append(ax_deep)    

        moc = MOC.sel(transport_reg=transport_reg)
        for ax in [ax_surf, ax_deep]:
            cs = ax.contour(
                moc.lat_aux_grid, moc.moc_z*1e-2, moc,
                colors='k',
                linewidths=1,
                levels=cnlevels,
            )

            cs.levels = [f'{val:g}' for val in cs.levels]
            fmt = '%r'
            lb = plt.clabel(cs, fontsize=8, inline=True, fmt=fmt)

            mesh = ax.contourf(
                moc.lat_aux_grid, moc.moc_z*1e-2, moc,
                levels=cnlevels,
                cmap='PRGn',
                extend='both',
                zorder=-1000,
            )

        ax_surf.set_ylim([1000., 0.])
        ax_surf.set_yticks(np.arange(0, 1000, 200))
        ax_surf.set_xticklabels([])
        ax_surf.xaxis.set_ticks_position('top')
        ax_surf.set_xticks(np.arange(-90, 110, 30))

        ax_deep.set_ylim([5000., 1000.])
        ax_deep.xaxis.set_ticks_position('bottom')
        ax_deep.set_xticks(np.arange(-90, 110, 30))

        if n == 0:
            ax_deep.set_ylabel('Depth [m]')
            ax_deep.yaxis.set_label_coords(-0.18, 1.05)
        else:
            ax_surf.set_yticklabels('')
            ax_deep.set_yticklabels('')

        ax_deep.set_xlabel('Latitude [°N]')

        ax_surf.set_title('MOC', loc='left')
        ax_surf.set_title(transport_reg, loc='center')    
        ax_surf.set_title('Sv', loc='right')                            
        ax_deep.set_xticklabels('')

    utils.label_plots(fig, [ax for ax in axs_surf] , xoff=-0.02, yoff=0.06)    

    
def tick_format(n: float, digits=None):
    if digits is None:
        digits = 5
    s = f'%.{digits}f' % n
    s += '.' if '.' not in s else ''
    if s[-1] == '0':
        return tick_format(n, digits=digits-1)
    else:
        return s[:-1] if s[-1] == '.' else s