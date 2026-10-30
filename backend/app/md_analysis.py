#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
FILES = {"rmsd":("rmsd.xvg","RMSD (nm)"),"rmsf":("rmsf.xvg","RMSF (nm)"),
         "gyrate":("gyrate.xvg","Radius of gyration (nm)"),"hbonds":("hbnum.xvg","Hydrogen bonds (count)")}
def read_xvg(path):
    pts=[]
    for line in path.read_text(errors="replace").splitlines():
        line=line.strip()
        if not line or line[0] in "#@&": continue
        a=line.split()
        try:
            x,y=float(a[0]),float(a[1])
            if math.isfinite(x) and math.isfinite(y): pts.append((x,y))
        except (ValueError,IndexError): pass
    return pts
def svg_plot(points,title,ylabel,path):
    w,h=900,420; left,right,top,bottom=85,25,55,65
    xs=[p[0] for p in points]; ys=[p[1] for p in points]
    xmin,xmax=min(xs),max(xs); ymin,ymax=min(ys),max(ys)
    if xmax==xmin: xmax=xmin+1
    if ymax==ymin: ymax=ymin+1
    pad=(ymax-ymin)*.08; ymin-=pad; ymax+=pad
    X=lambda x:left+(x-xmin)/(xmax-xmin)*(w-left-right)
    Y=lambda y:h-bottom-(y-ymin)/(ymax-ymin)*(h-top-bottom)
    coords=" ".join(f"{X(x):.1f},{Y(y):.1f}" for x,y in points)
    esc=lambda x:str(x).replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")
    svg=(f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">'
         f'<rect width="100%" height="100%" fill="white"/>'
         f'<text x="{w/2}" y="30" text-anchor="middle" font-family="Arial" font-size="20">{esc(title)}</text>'
         f'<line x1="{left}" y1="{top}" x2="{left}" y2="{h-bottom}" stroke="black"/>'
         f'<line x1="{left}" y1="{h-bottom}" x2="{w-right}" y2="{h-bottom}" stroke="black"/>'
         f'<polyline points="{coords}" fill="none" stroke="#1769aa" stroke-width="2"/>'
         f'<text x="{w/2}" y="{h-15}" text-anchor="middle" font-family="Arial" font-size="14">X value (as reported by GROMACS)</text>'
         f'<text x="20" y="{h/2}" transform="rotate(-90 20 {h/2})" text-anchor="middle" font-family="Arial" font-size="14">{esc(ylabel)}</text>'
         f'<text x="{left}" y="{h-bottom+22}" font-family="Arial" font-size="11">{xmin:.3g}</text>'
         f'<text x="{w-right}" y="{h-bottom+22}" text-anchor="end" font-family="Arial" font-size="11">{xmax:.3g}</text></svg>')
    path.write_text(svg)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--dir",default="."); d=Path(ap.parse_args().dir)
    out={"files":{},"note":"Descriptive XVG summaries only; these do not validate MD convergence."}
    for key,(filename,label) in FILES.items():
        p=d/filename
        if not p.exists(): out["files"][key]={"status":"missing"}; continue
        pts=read_xvg(p)
        if not pts: out["files"][key]={"status":"no_numeric_data"}; continue
        ys=[y for _,y in pts]
        out["files"][key]={"status":"ok","n_points":len(pts),"x_start":pts[0][0],"x_end":pts[-1][0],
          "min":min(ys),"max":max(ys),"mean":sum(ys)/len(ys),"final":ys[-1],"source":filename,"plot":f"{key}.svg"}
        svg_plot(pts,key.upper(),label,d/f"{key}.svg")
    (d/"md_analysis_summary.json").write_text(json.dumps(out,indent=2))
    with (d/"md_analysis_summary.csv").open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["metric","status","n_points","min","max","mean","final","source"])
        for k,v in out["files"].items(): w.writerow([k,v.get("status"),v.get("n_points"),v.get("min"),v.get("max"),v.get("mean"),v.get("final"),v.get("source")])
    print("Wrote MD analysis CSV/JSON and SVG plots.")
if __name__=="__main__": main()
