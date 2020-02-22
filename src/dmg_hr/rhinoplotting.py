import json
import rhinoscriptsyntax as rs

from compas.utilities.colors import i_to_rgb

rs.DeleteObjects(rs.ObjectsByLayer('Default'))
filepath = '/Users/time/Documents/UW/04_code/dmg_helmholtz/data/panel_geometry.json'

with open(filepath, 'r') as fp:
    data = json.load(fp)

rs.AddTextDot('S', data['src'])

hr_pts = data['hr_pts']
nr = data['nr']
nl = data['nl']
br = data['br']
bl = data['bl']
for i, hr_pt in enumerate(hr_pts):
    r = nr[i]
    npt = [hr_pt[0], hr_pt[1], hr_pt[2] - nl[i]]
    # bpt = [hr_pt[0], hr_pt[1], hr_pt[2] - nl = bl]
    rs.AddCylinder(hr_pt, -nl[i], r, cap=False)
    rs.AddCylinder(npt, -bl[i], br[i], cap=False)

recs = data['recs']
ps = data['pressures']
minp = min(ps)
maxp = max(ps)
for i, rec in enumerate(recs):
    pt = rs.AddPoint(rec)
    p = ps[i]
    i = (((p - minp) * (1 - 0)) / (maxp - minp)) + 0
    color = i_to_rgb(i)
    rs.ObjectColor(pt, color)
