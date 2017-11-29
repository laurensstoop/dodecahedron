#!/usr/bin/env python3

import numpy as np
#import math
import h5py

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
from matplotlib.path import Path
from matplotlib.transforms import Affine2D
import sys

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import matplotlib.colors as colors


# The scripts 
sys.path.append('/usr/people/stoop/Documents/thesis/scripts') 
from definitions_general import ReadData, PlotMap
from definitions_windenergy import WindAtHeight, NormalizedWindPowerCurve


class dummy_with(object):
    def __init__(self):
        pass
    def __enter__(self):
        pass
    def __exit__(self, t, value, tb):
        pass
    def savefig(self):
        pass

def dodecahedron_cornerpoints():
    a = np.arccos((np.sqrt(5)-1)/(2*np.sqrt(3)*np.sin(np.pi/5)))
    b = a - 2*np.arcsin((np.sqrt(5)-1)/(2*np.sqrt(3)))
    c = np.pi/2 - a - b

    a *= 180/np.pi
    b *= 180/np.pi
    c *= 180/np.pi

    coordinates = [{'lat': np.zeros((6,), dtype=np.float64) + a,
                    'lon': np.arange(36.0, 400.0, 72.0),
                    'lat_c': 90.0, 'lon_c': 0.0},
                   {'lat': np.asarray([a, a, b, -b, b, a]),
                    'lon': np.asarray([180.0, 252.0, 252.0, 216.0, 180.0, 180.0]),
                    'lat_c': c, 'lon_c': 216.0},
                   {'lat': np.asarray([a, a, b, -b, b, a]),
                    'lon': np.asarray([252.0, 324.0, 324.0, 288.0, 252.0, 252.0]),
                    'lat_c': c, 'lon_c': 288.0},
                   {'lat': np.asarray([a, a, b, -b, b, a]),
                    'lon': np.asarray([324.0, 36.0, 36.0, 0.0, 324.0, 324.0]),
                    'lat_c': c, 'lon_c': 0.0},
                   {'lat': np.asarray([a, a, b, -b, b, a]),
                    'lon': np.asarray([36.0, 108.0, 108.0, 72.0, 36.0, 36.0]),
                    'lat_c': c, 'lon_c': 72.0},
                   {'lat': np.asarray([a, a, b, -b, b, a]),
                    'lon': np.asarray([108.0, 180.0, 180.0, 144.0, 108.0, 108.0]),
                    'lat_c': c, 'lon_c': 144.0},
                   {'lat': np.asarray([b, -b, -a, -a, -b, b]),
                    'lon': np.asarray([252.0, 288.0, 288.0, 216.0, 216.0, 252.0]),
                    'lat_c': -c, 'lon_c': 252.0},
                   {'lat': np.asarray([b, -b, -a, -a, -b, b]),
                    'lon': np.asarray([324.0, 0.0, 0.0, 288.0, 288.0, 324.0]),
                    'lat_c': -c, 'lon_c': 324.0},
                   {'lat': np.asarray([b, -b, -a, -a, -b, b]),
                    'lon': np.asarray([36.0, 72.0, 72.0, 0.0, 0.0, 36.0]),
                    'lat_c': -c, 'lon_c': 36.0},
                   {'lat': np.asarray([b, -b, -a, -a, -b, b]),
                    'lon': np.asarray([108.0, 144.0, 144.0, 72.0, 72.0, 108.0]),
                    'lat_c': -c, 'lon_c': 108.0},
                   {'lat': np.asarray([b, -b, -a, -a, -b, b]),
                    'lon': np.asarray([180.0, 216.0, 216.0, 144.0, 144.0, 180.0]),
                    'lat_c': -c, 'lon_c': 180.0},
                   {'lat': np.zeros((6,), dtype=np.float64) - a,
                    'lon': np.arange(0.0, 400.0, 72.0),
                    'lat_c': -90.0, 'lon_c': 0.0}]
    return coordinates

def dodecahedron_globe(species=None, data=None, format='svg', indices=None, bar=False):
    if species is None:
        fname = './img/dodecahedron_globe.pdf'
        base_fname = './img/dodecahedron_globe_{{0:02d}}.{0}'.format(format)
    else:
        fname = './img/dodecahedron_globe_{0}.pdf'.format(species)
        base_fname = './img/dodecahedron_globe_{0}_{{0:02d}}.{1}'.format(species, format)

    if indices is None:
        indices = range(12)
        pdf = PdfPages(fname)
    else:
        indices = [i-1 for i in indices]
        pdf = dummy_with()
    coordinates = dodecahedron_cornerpoints()

    with pdf:
        for idx in indices:
            d = coordinates[idx]
            print("figure {0}".format(idx+1))
            print("    center   : {0:.2f}, {1:.2f}".format(d['lon_c'], d['lat_c']))
            print("    longitude: [{0[0]:.2f}, {0[1]:.2f}, {0[2]:.2f}, {0[3]:.2f}, {0[4]:.2f}]".format(d['lon']))
            print("    latitude : [{0[0]:.2f}, {0[1]:.2f}, {0[2]:.2f}, {0[3]:.2f}, {0[4]:.2f}]".format(d['lat']))

            fig = plt.figure(0, figsize=(7.5, 7.5))
            m = Basemap(projection='gnom',
                        width=1e7, height=1e7,
                        lat_0=d['lat_c'],lon_0=d['lon_c'],
                        resolution='i')
            if species is None:
                m.drawmapboundary(fill_color='lightskyblue', zorder=1)
                m.fillcontinents(color='darkgreen',lake_color='lightskyblue', zorder=2)
                # m.drawparallels(np.arange(-80, 90, 10), zorder=3, color='white')
                # m.drawmeridians(np.arange(0.0, 360.0, 18.0), zorder=3, color='white')
                # m.drawrivers(linewidth=0.5, linestyle='solid', color='lightskyblue', zorder=3)
                # m.drawcoastlines(linewidth=0.25, color='w', zorder=3, linestyle='solid')
            else:
                m.drawcoastlines(linewidth=0.5, color='w', zorder=3, linestyle='solid')

            # plot data.
            if species is not None:
                print("    Plotting data for {0}".format(species))
                cmap = cm.get_cmap('BuGn')#name='jet' if species is 'NO2' else 'ozone_special')
                if species == 'O3':
                    cmap.set_under('b')
                    cmap.set_over('w')

                if species is 'NO2':
                    normalizer = colors.PowerNorm(vmin=data['range'][0],
                                                  vmax=data['range'][1],
                                                  gamma=1/2,
                                                  clip=True)
                else:
                    normalizer = colors.Normalize(vmin=data['range'][0],
                                                  vmax=data['range'][1],
                                                  clip=True)


                lons, lats = np.meshgrid(data['lon'], data['lat'])
                
                if species is 'rgb':                
                    # The following pcolormesh code does work when I want to plot the data to normal full world map
                    # projections like mollweide or robinson projections, but unfortunately it fails here.
                    #cs = m.pcolormesh(lons, lats, data['var'][:,:,0], latlon=True, zorder=1,
                    #            color=np.array([data['var'][:, :, 0].ravel(), data['var'][:, :, 1].ravel(), data['var'][:, :, 2].ravel()]).transpose()).set_array(None)
                    #cs = m.pcolor(lons, lats, data['var'][:,:,0], latlon=True, zorder=1,
                    #                color=np.array([data['var'][:, :, 0].ravel(), data['var'][:, :, 1].ravel(), data['var'][:, :, 2].ravel()]).transpose())

                    # Therefore I utilised the basemap transform function in combination with imshow,
                    # which since basemap version 1.1.0 does not contain a bug anymore.                    
                    m.imshow(np.dstack([m.transform_scalar(data['var'][:, :, c], data['lond'], data['latd'], 
                                                           data['plon'], data['plat']) for c in range(3)]))
                                    
                else:
                    cs = m.pcolor(lons, lats, data['var'], cmap=cmap, latlon=True, edgecolors='face',
                                  zorder=1, norm=normalizer)
                if bar:
                    ticks = [1e15,5e15,1e16,1.5e16]
                    cbar = m.colorbar(cs, location='bottom', pad="5%", extend="both",
                                      ticks=ticks if species == 'NO2' else None)
                    labeltxt = ''
                    if species == 'O3':
                        labeltxt = "Total ozone column [DU]" 
                    elif species == 'NO2':
                        labeltxt = 'NO$_2$ Concentration [molecules/cm$^2$]'
                    elif species == 'EC-Earth':
                        labeltxt = 'E$_{out}$ [GWh]'
                        
                    cbar.set_label(labeltxt)
                    
            print("    Adding pentagon")
            x, y = m(d['lon'], d['lat'])
            xy = np.concatenate((x, y)).reshape((2,6)).T

            poly = ptch.Polygon(xy, closed=True, fc='none', ec='r', zorder=4)
            # poly = Path(xy, closed=True)
            ax = plt.gca()
            ax.add_patch(poly)

            poly_k = ptch.Polygon(xy, closed=True, fc='none', ec='k', zorder=5)
            ax.add_patch(poly_k)

            if idx == 0:
                xx, yy = m(0.0, 90.0)
                plt.plot(xx, yy, 'wx')
            elif idx == 11:
                xx, yy = m(0.0, -90.0)
                plt.plot(xx, yy, 'wx')

            if species is None:
                plt.title('figure {0}'.format(idx+1))
            else:
                plt.title('{0} part {1}'.format(species, idx+1))

            print("    Saving figure to {0}".format(base_fname.format(idx+1)))
            fig.savefig(base_fname.format(idx+1), dpi=600)

            print("    Adding figure to pdf output")
            pdf.savefig()
            plt.close(fig)

    if species is None:
        print("Done!")
    else:
        print("Done with {0}!".format(species))


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='Plot on a dodecahedron')
    parser.add_argument("--omi_rgb", action='store_true', help="Plot omi rgb")
    parser.add_argument("--ozone", action='store_true', help="Plot ozone")
    parser.add_argument("--ecearth", action='store_true', help="Plot EC-Earth")
    parser.add_argument("--no2", action='store_true', help="Plot NO2")
    parser.add_argument("--no2-full", action='store_true', help="Plot NO2 (full res).")
    parser.add_argument("--blank", action='store_true', help="Plot blank earth")
    parser.add_argument("--index", nargs='+', type=int,
                        choices=(1,2,3,4,5,6,7,8,9,10,11,12),
                        help="Specify index of image to plot")
    parser.add_argument('--format', choices=('png', 'svg', 'pdf'), default='png', help="output format")
    parser.add_argument('--bar', action='store_true', help="add color bar")

    args = parser.parse_args()

    # Omi rgb data is on a rectanguler lat,lon hd grid,
    # with a third dimension of size 3 for r, g and b values ranged from 0 to 1.
    if args.omi_rgb:
    
        print("Preparing OMI rgb data")
        data = {}
        with h5py.File('rgb_img_year_1080x1920.h5', 'r') as ref:
            var = np.asarray(ref['rgb_img'][:], dtype=np.float64)
        print("    Shape of data: [{0[0]}, {0[1]}]".format(var.shape))
        
        data['var'] = var        
        lat0, lat1 = -90.0, 90.0
        lon0, lon1 = -180.0, 180.0
        
        data['lon'] = np.linspace(lon0, lon1, num=var.shape[1])
        data['lat'] = np.linspace(lat0, lat1, num=var.shape[0])
        data['plon'], data['plat'] = var.shape[1], var.shape[0]
        dlon, dlat = (lon1 - lon0) / data['plon'], (lat1 - lon0) / data['plat']
        data['lond'] = np.linspace(lon0 + dlon, lon1 - dlon, data['plon'])
        data['latd'] = np.linspace(lat0 + dlat, lat1 - dlat, data['plat'])
        data['range'] = [0.0, 1.0]
        dodecahedron_globe(species='rgb', data=data, format=args.format, indices=args.index, bar=args.bar)
	
    if args.ozone:
        print("Preparing Ozone data")
        data = {}

        # The file 'o3col2006102112.h5' is converted from the original
        # 'o3col2006102112.hdf' using h4toh5. The original fil is available
        # on TEMIS (http://temis.nl/protocols/o3field/o3field_msr2.php?Year=2006&Month=10&Day=13)
        with h5py.File('o3col2006102112.h5', 'r') as ref:
            var = np.asarray(ref['O3_column'][:], dtype=np.float32)

        print("    Shape of data: [{0[0]}, {0[1]}]".format(var.shape))
        data['var'] = var
        data['lon'] = np.linspace(-180.0, 180.0, num=var.shape[1])
        data['lat'] = np.linspace(-90.0, 90.0, num=var.shape[0])
        data['range'] = [100.0, 450.0]

        dodecahedron_globe(species='O3', data=data, format=args.format, indices=args.index, bar=args.bar)

    if args.no2_full or args.no2:
        data = {}

        # This data file was provided by Pepijn veefkind
        with h5py.File('domino_monthly_climatology_2005_2016_0.125x0.125.h5', 'r') as ref:
            if args.no2_full:
                print("Preparing full res NO2 data")
                var = np.nanmean(ref['NO2'][:,:,:], axis=0)
            else:
                print("Preparing sliced NO2 data for reduction by factor 2")
                var0 = np.nanmean(ref['NO2'][:,0::2,0::2], axis=0)
                var1 = np.nanmean(ref['NO2'][:,0::2,1::2], axis=0)
                var2 = np.nanmean(ref['NO2'][:,1::2,0::2], axis=0)
                var3 = np.nanmean(ref['NO2'][:,1::2,1::2], axis=0)
                var = np.nanmean(np.asarray([var0, var1, var2, var3]), axis=0)

                del var0, var1, var2, var3

        print("    Shape of data: [{0[0]}, {0[1]}]".format(var.shape))
        data['var'] = var
        data['lon'] = np.linspace(-180.0, 180.0, num=var.shape[1])
        data['lat'] = np.linspace(-90.0, 90.0, num=var.shape[0])
        data['range'] = [5.0e14, 1.5e16]

        dodecahedron_globe(species='NO2', data=data, format=args.format, indices=args.index, bar=args.bar)


    if args.ecearth:
        print("Preparing EC-Earth data")
        data = {}
        
        # Get a datahandler
        import xarray as xr
                
        # Open the dataset        
        data_file = xr.open_dataset('sfcwind_d_ECEarth_PD_s01r00_2035_coord.nc')
        
        # Calculate the wind at height 80
        data_file['wind80'] = WindAtHeight(data_file,80)
        
        # Calculate the normalized wind power curve
        data_file['wind80powercurve'] = NormalizedWindPowerCurve(data_file)
        
        # Calculate the mean
        data_mean_t = data_file.mean(dim='time', keep_attrs=True)
        
        # Wind energy density (W m**-2) As (Hueging et al 2012)
        data_mean_t['wed80'] = 0.5*1.225*data_mean_t['wind80']**3
        
        # Extractable wind energy for 3MW windmill per gridcell
        data_file['windenergy80'] = 3e6 * data_file['wind80powercurve'] / 1e9 # GW
        
        # Calculate the sum over the whole year
        data_sum_t = data_file.sum(dim='time')
        
        
        # Calculate the mean
        data_mean_t = data_file.mean(dim='time', keep_attrs=True)
            
        # Old stuff
        #print("    Shape of data: [{0[0]}, {0[1]}]".format(var.shape))
        #data['var'] = data_mean_t['wed80']
        data['var'] = data_sum_t['windenergy80']
        
        data['lon'] = np.linspace(-180.0, 180.0, num=320)
        data['lat'] = data_mean_t['lat']#np.linspace(-90.0, 90.0, num=var.shape[0])
        data['range'] = [0.0, 1.0]

        dodecahedron_globe(species='EC-Earth', data=data, format=args.format, indices=args.index, bar=args.bar)
    
    
    if args.blank:
        dodecahedron_globe(format=args.format, indices=args.index)
