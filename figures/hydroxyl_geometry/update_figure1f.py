"""Regenerate only Figure 1F's illustrative reference-width guides.

The scientific source generator and CSV files remain unchanged. PNG outputs only.
"""
from pathlib import Path
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parent
WORK=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager
from matplotlib.transforms import Bbox
import numpy as np
from PIL import Image,ImageChops
import msfig_common as common

for font in Path('/mnt/c/Windows/Fonts').glob('arial*.ttf'):
    font_manager.fontManager.addfont(str(font))
common.TAB=str(WORK)
source_path=WORK/'msfig1.py'
source=source_path.read_text(encoding='utf-8')
# RDKit is unnecessary for a panel-F-only revision: panels A–D are retained
# directly from the current document artwork, while E is also preserved.
for line in ['from rdkit import Chem\n','from rdkit.Chem import AllChem\n',
             'from rdkit.Chem.Draw import rdMolDraw2D\n']:
    source=source.replace(line,'')
start=source.index('    for k, s in enumerate(SYS):')
end=source.index('    # ---- E:',start)
source=source[:start]+source[end:]
old="""    axF.axhline(8.6, color='black', ls=(0, (4, 2)), lw=0.8)
    axF.text(-0.42, 8.9, 'reach needed to bridge\\nboth ends (8.6 Å)',
             ha='left', va='bottom', fontsize=6.3, linespacing=1.15)"""
new="""    # Illustrative fixed-width guides, not universal bridging cutoffs.
    axF.axhline(8.62, color='black', ls=(0, (4, 2)), lw=0.8)
    axF.axhline(7.22, color='#666666', ls=(0, (1, 2)), lw=0.8)
    axF.text(-0.42, 11.8, 'Reference W = 14.22 Å',
             ha='left', va='bottom', fontsize=6.3)
    axF.text(-0.42, 8.9, 'W − 2a = 8.62 Å (Eq. 3)',
             ha='left', va='bottom', fontsize=6.3)
    axF.text(-0.42, 7.45, 'W − 2$d_0$ = 7.22 Å (Eq. 4)',
             ha='left', va='bottom', fontsize=6.3, color='#666666')"""
assert source.count(old)==1
records={}
panel_bounds={}
print_bbox=None

def geometry(fig):
    return [{'position':list(ax.get_position().bounds),
             'patches':[[p.get_x(),p.get_y(),p.get_width(),p.get_height()]
                        for p in ax.patches if hasattr(p,'get_height')],
             'images':[np.asarray(i.get_array()).tolist() for i in ax.images],
             'xlim':list(ax.get_xlim()),'ylim':list(ax.get_ylim())} for ax in fig.axes]

for version,code in [('baseline',source),('study',source.replace(old,new))]:
    scope={'__name__':'study_figure1_generator','__file__':str(source_path)}
    exec(compile(code,str(source_path),'exec'),scope)
    def save(fig,ignored_name,version=version):
        global print_bbox
        records[version]=geometry(fig)
        fig.set_dpi(600)
        fig.canvas.draw()
        if print_bbox is None:
            lower=fig.get_tightbbox(fig.canvas.get_renderer()).padded(0.02)
            print_bbox=Bbox.from_bounds(lower.x0,lower.y0,4112/600,2100/600)
        bounds=fig.axes[-1].get_window_extent(fig.canvas.get_renderer())
        panel_bounds[version]=[bounds.x0-print_bbox.x0*600,print_bbox.y1*600-bounds.y1,
                               bounds.x1-print_bbox.x0*600,print_bbox.y1*600-bounds.y0]
        fig.savefig(WORK/f'Figure1F_{version}_canvas.png',bbox_inches=print_bbox,pad_inches=0)
    scope['save']=save
    scope['main']()
assert records['baseline']==records['study'],'Geometry or underlying plot data changed'
original=Image.open(WORK/'Figure1_original.png').convert('RGBA')
baseline=Image.open(WORK/'Figure1F_baseline_canvas.png').convert('RGBA')
canvas=Image.open(WORK/'Figure1F_study_canvas.png').convert('RGBA')
def color_bbox(im,color=(0,114,178,255)):
    points=np.argwhere((np.asarray(im)==np.array(color)).all(axis=2))
    return (int(points[:,1].min()),int(points[:,0].min()),int(points[:,1].max()+1),int(points[:,0].max()+1))
ob,bb=color_bbox(original),color_bbox(baseline)
assert ob[2]-ob[0]==bb[2]-bb[0] and ob[3]-ob[1]==bb[3]-bb[1],(ob,bb)
dx,dy=ob[0]-bb[0],ob[1]-bb[1]
box=tuple(round(v) for v in panel_bounds['baseline'])
# Keep the original black axes intact and update only their interior. Every
# baseline pixel in that region must first reproduce the original exactly.
box=(box[0]+2,box[1]+2,box[2]-2,box[3]-2)
target_box=tuple(v+(dx if i%2==0 else dy) for i,v in enumerate(box))
baseline_bbox=ImageChops.difference(original.crop(target_box).convert('RGB'),baseline.crop(box).convert('RGB')).getbbox()
assert baseline_bbox is None,('Original F plot did not regenerate identically',baseline_bbox,target_box,box)
updated=original.copy()
updated.paste(canvas.crop(box),target_box)
updated.save(WORK/'Figure1.png')
diff_bbox=ImageChops.difference(original.convert('RGB'),updated.convert('RGB')).getbbox()
assert diff_bbox is not None and diff_bbox[0]>2000 and diff_bbox[1]>900, diff_bbox
# The region through panels A–E must remain pixel-identical; all differences
# are constrained to panel F's former annotation/guide area.
assert ImageChops.difference(original.crop((0,0,2000,2100)).convert('RGB'),updated.crop((0,0,2000,2100)).convert('RGB')).getbbox() is None
assert ImageChops.difference(original.crop((0,0,4112,900)).convert('RGB'),updated.crop((0,0,4112,900)).convert('RGB')).getbbox() is None
data={
 'source_generator':source_path.name,
 'source_generator_sha256':hashlib.sha256(source_path.read_bytes()).hexdigest(),
 'source_image_sha256':hashlib.sha256((WORK/'Figure1_original.png').read_bytes()).hexdigest(),
 'output':'Figure1.png',
 'output_sha256':hashlib.sha256((WORK/'Figure1.png').read_bytes()).hexdigest(),
 'dimensions':updated.size,
 'regenerated_panel_F_baseline_pixel_identical':True,
 'panels_A_to_E_pixel_identical':True,
 'underlying_bar_data_and_geometry_unchanged':True,
 'changed_pixel_bbox':diff_bbox,
 'panel_F_interior_box':target_box,'canvas_to_source_translation':[dx,dy],
 'reference_width_A':14.22,'exact_source_reference_width_A':14.216,
 'allowance_A':2.8,'contact_cutoff_A':3.5,
 'margin_zero_guide_A':8.62,'contact_bound_guide_A':7.22,
 'source_csv':'s4_an_geometry.csv',
 'reference_source_row':'Anthraflavic acid (2,6-diOH)',
 'reference_source_field':'pocket_GLU353_HIS524_mean_A',
 'plot_interpretation':'Illustrative guides at the rounded first-2,6-run mean width; not universal cutoffs.',
 'pdf_generated':False,
}
(WORK/'figure1f_validation.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(data,ensure_ascii=False,indent=2))
