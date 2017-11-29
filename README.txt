Dodecahedron Earth

How to make your own globe on a dodecahedron, including earth observation data.

Tools used:

* Python 3, with numpy, matplotlib, Basemap and h5py.
* InkScape (or Adobe Illustrator)

All initial steps use the dodecahedron_earth.py script.

The template 'dodecahedron_glob_template.svg' was created from the 'blank'
output.

    python3 dodecahedron_earth.py --blank --format svg

This will generate sgv files without earth observation data. These can be used
in the template to ensure everything is positioned correctly. The squares in the
template are the size of the png files that come out of the other calls. These
can be used to position those png files exactly where needed. It is not possible
(in my experience) to use svg or pdf output from the dodecahedron_earth.py for
earth observation data directly. The resulting files are simply too large. In
the template all land outlines were removed, only the square outlining the png
and the pentagon to cut out the required area are left. The numbers 1-12
correspond with the parts, the +36 and -36 figures indicate the rotation in
degrees that needs to be applied to that figure (negative is clockwise).


The next step is to create image files with earth observation data.

    python3 dodecahedron_earth.py --ozone --format png

This will generate a sequence of png files that can be pasted in (a copy of!)
the template. The other data for which the script is prepared is '--no2' or
'--no2-full' (both use the same data, but one downsamples the data to reduce
plotting time).

If you want a color scale then call the script again (after you've copied the
png files into your copy of the template):

    python3 dodecahedron_earth.py --ozone --format png --bar --index 1

This will re-create the first output file, but with a color bar.


The code in the script is simple enough (I think) to add your own data. Do not
try to further optimize the plotting code with pcolormesh(), that won't work.

When EC-Earth(like) data is used (so --ecearth instead of --ozone), just copy the template into a folder containing the 12 figures and 1 bar figure. Open the template then in inktscape and save as pdf. See folders in ./img/ for more detail. 


The square outlines allow for exact positioning, and the pentagons can be used
to clip the png files to show only what is needed for the papercraft.
In InkScape the template can be used with the following steps:

1. Make a copy of the template.
2. Lock all layers (show the layers from the Layers menu, using the 'Layers...'
   menu item). Unlock and select the Pentagon layer. You can hide the
   "Plakstroken" layer (the grey tabs for gluing).
3. Import your png file. Move the image to the correct position. InkScape will
   snap to the corners of the squares (and the pentagons, make sure you are
   locked correctly).
4. Lower the png image within the pentagon layer, such that the pentagon is on
   top of the png image.
5. Select the png image and then (while holding the shift key) select the
   pentagon. Both should now be selected.
6. Choose 'Set' from the 'Clip' sub-menu in the 'Object' menu.
7. Repeat steps 3-6 for all 11 other png files. Make sure to rotate as
   indicated, using the 'Transform...' menu item in the 'Object' menu.
8. Add extra items needed for the hand-out (copyright, description color bar,
   ...)
9. Hide the Squares and Numbers layers and save a copy of the file as pdf.
   This file can be printed.

The dodecahedron_earth.py script also produces pdf files with full resolution
output on 12 pages. These can be used to make larger versions of the
dodecahedrons.

Maarten Sneep (maarten.sneep@knmi.nl), 2017-10-04.
Laurens Stoop (laurensstoop@protonmail.com), 2017-11-29.
