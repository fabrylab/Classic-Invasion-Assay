from clickpoints import DataFile
import glob
import os
import peewee
import sys


# This programm is used to generate .cdb clickpoints databases. Image z-stacks with BF and Fluo images are required.
# Der Modus (BF/Fluo) dient als Layer und die eigentliche Ebene (layer) vom z-Stack wird als Index (frame number)
# in die Datenbank einsortiert. So können die Daten für den Invasions-Assay mit Annalena evaluiert werden in Clickpoints.

src_path = r'H:\Experiment_data\B01_AnWi_Invasion_2020-01-31\U87-1\pos001'

if len(sys.argv) == 2:
    src_path = sys.argv[1]
    print("using path: %s" % src_path)


channels = ['modeBF','modeFluo5']


# get files
bf_images = glob.glob(os.path.join(src_path, '*' + channels[0] + '*'))
fluo_images = glob.glob(os.path.join(src_path, '*' + channels[1] + '*'))

db_name = os.path.join(src_path, 'sorted.cdb')
if os.path.exists(db_name):
    os.remove(db_name)
    db = DataFile(db_name,'w')
else:
    db = DataFile(db_name,'w')

c0 = db.setLayer(channels[0])
c1 = db.setLayer(channels[1])


for im in bf_images:
    path,file = os.path.split(im)


    p = db.setPath('.')
    try:
        db.setImage(file,p,layer=c0)
    except peewee.IntegrityError:
        pass

for id,im in enumerate(fluo_images):
    path, file = os.path.split(im)

    p = db.setPath('.')
    try:
        db.setImage(file, p, layer=c1, sort_index=id)
    except peewee.IntegrityError:
        pass


db.setMarkerType(name='cell_in_focus', color= '#ff0000', mode=db.TYPE_Track)
db.setMarkerType(name='not_a_cell', color= '#e2ff00', mode=db.TYPE_Track)

db.db.close()