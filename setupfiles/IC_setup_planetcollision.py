## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
import IC_denstable as denstable
from scipy.integrate import simps, odeint  # Add simps to imports
from scipy.interpolate import interp1d
#from createplanet_2025_b import get_inputs, run_model
from createplanet_2025_b import run_or_load
import re

class setup_planetcollision(object):
    dICdensRsmooth = 0.1
    dICdensprofile = 8
    dICdensdir = 4
    dICdensR = 1.0 
    dICdensinner = 1.0 # calculated depending on mass
    dICdensouter = 1.0
    inflow=0
    rhoit=0
    npart=0
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    ns=64
    time=0
    grav=1
    cosmo=0
    molweight=1.0
    gamma=1.6666666666667
    periodic=0
    deltastep=0.1
    freqout=1
    nsteps=3000
    dmsolunit=1.0
    dkpcunit=1.0
    adi=2
    shape=2
    mass=[]
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]
    rho=[]
    u=[]
    h=[]
    P=[]
    Bx=[]
    By=[]
    Bz=[]
    Mat1=[]
    Mat2=[]
    Mat3=[]
    Mat4=[]
    Mat5=[]
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    M=1.0
    R=1.0
    Rtab = 1.0
    eps = 0.01
    materialfile=""
   
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='glassreadyfile',target=None,impactor=None,impang=0.73,vinf=0.0,impgeom=2):
         if target == None or impactor == None or target == "None" or impactor == "None":
             print("Specify target and impactor")
             exit
         # We want stuff in units of Earth mass and Earth Radius
         M_earth = 5.972*10**24 #kg
         R_earth = 6378*1000 # m
         M_sol = 1.989*10**30 #kg
         R_kpc = 3.086*10**19 #m
         #we set M_earth and R_earth as code units.
         #1 Code unit =  1 Msol unit * M_earth/M_sol
         self.dmsolunit = M_earth/M_sol
         self.dkpcunit = R_earth/R_kpc       
         self.nsteps = round(48.0/(806.727*self.deltastep/(60*60))) #48hrs / hr per dump
         Bzero = 0.0;
         tgdatatar,tddata,tsdata,data_header,time=tip.readtipsy(target);
         ngastar = data_header[2]
         Xtar,labelstar=tip.readalltipsyaux(target)
         tar_matfile=self.getmaterialfilename(target)
         tgdataimp,tddata,tsdata,data_header,time=tip.readtipsy(impactor);
         ngasimp = data_header[2]
         Ximp,labelsimp=tip.readalltipsyaux(impactor)
         imp_matfile=self.getmaterialfilename(impactor)

         Rtar=np.max(np.sqrt(tgdatatar[:,1]**2+tgdatatar[:,2]**2+tgdatatar[:,3]**2))
         Rimp=np.max(np.sqrt(tgdataimp[:,1]**2+tgdataimp[:,2]**2+tgdataimp[:,3]**2))
         Mimp=np.sum(tgdataimp[:,0])
         Mtar=np.sum(tgdatatar[:,0])
         dist=(Rtar+Rimp)*1.001
         tgdataimp[:,1]=tgdataimp[:,1]+dist
         print(dist,impgeom,ngasimp,ngastar)
         impgeom=int(impgeom)
         impang=float(impang)
         vinf=float(vinf)
         tgdataimp[:,impgeom]=tgdataimp[:,impgeom]+dist*impang
         G=1.0
         vesc=np.sqrt((2*G*(Mtar+Mimp))/(Rtar+Rimp))
         vimp=np.sqrt(vesc*vesc+vinf*vinf)

         vximp =  -(Mtar/(Mtar+Mimp)) * vimp
         vxtar =  (Mimp/(Mtar+Mimp)) * vimp

         cmx = (np.sum(tgdataimp[:,0]*tgdataimp[:,1])+np.sum(tgdatatar[:,0]*tgdatatar[:,1]))/(Mtar+Mimp)
         cmy = (np.sum(tgdataimp[:,0]*tgdataimp[:,2])+np.sum(tgdatatar[:,0]*tgdatatar[:,2]))/(Mtar+Mimp)
         cmz = (np.sum(tgdataimp[:,0]*tgdataimp[:,3])+np.sum(tgdatatar[:,0]*tgdatatar[:,3]))/(Mtar+Mimp)

         tgdataimp[:,1]=tgdataimp[:,1]-cmx
         tgdataimp[:,2]=tgdataimp[:,2]-cmy
         tgdataimp[:,3]=tgdataimp[:,3]-cmz

         tgdatatar[:,1]=tgdatatar[:,1]-cmx
         tgdatatar[:,2]=tgdatatar[:,2]-cmy
         tgdatatar[:,3]=tgdatatar[:,3]-cmz

         tgdataimp[:,4]=tgdataimp[:,4]+vximp
         tgdatatar[:,4]=tgdatatar[:,4]+vxtar
         
         matfile="collision_material"
         self.materialfile=matfile
         tarmatidx,impmatidx=self.combine_material_files_dedup(tar_matfile,imp_matfile,matfile)
         print("TEST",tarmatidx,impmatidx)
         # Set materials in correct way:
         orig_tar = Xtar.copy()
         orig_imp = Ximp.copy()
         for i in range(len(tarmatidx)):
             orgi=labelstar.index(f"Material{i+1}")
             newi=labelstar.index(f"Material{tarmatidx[i]+1}")
             Xtar[:,newi]=orig_tar[:,orgi]

         for i in range(len(impmatidx)):
             orgi=labelsimp.index(f"Material{i+1}")
             newi=labelsimp.index(f"Material{impmatidx[i]+1}")
             Ximp[:,newi]=orig_imp[:,orgi]


         # Concatenate impactor and target data
         self.Xg = np.concatenate((Xtar, Ximp));
         self.labels = labelstar;
         tgdata = np.concatenate((tgdatatar, tgdataimp))
         self.mass=tgdata[:,0]
         self.x=tgdata[:,1]
         self.y=tgdata[:,2]
         self.z=tgdata[:,3]
         self.vx=tgdata[:,4]
         self.vy=tgdata[:,5]
         self.vz=tgdata[:,6]
         self.rho=tgdata[:,7]
         self.u=tgdata[:,8]
         #self.h=tgdata[:,9]
         self.npart=len(self.x)
         self.ngas=self.npart
         print('npart = ',self.npart,' particle mass = ',self.mass[1])
         self.h = [self.eps]*self.npart
         for i in range(self.npart):
             self.Bx.append(0.)
             self.By.append(0.)
             self.Bz.append(Bzero)

    def getmaterialfilename(self,filename):
            filename = re.sub(r"\.\d+$", "", filename)
            paramfile = f"{filename}.param"
            achMaterialFile = None

            with open(paramfile, "r") as f:
                for line in f:
                    if line.strip().startswith("achMaterialFile"):
                        # Split at '=' and take the right-hand side, strip whitespace
                        achMaterialFile = line.split("=", 1)[1].strip()
                        break

            print("achMaterialFile =", achMaterialFile)
            return achMaterialFile


    def combine_material_files(self,matfile1, matfile2, outfile):
        """
        Combine two material files into one.
        - Materials in matfile1 keep their indices (0,1,...).
        - Materials in matfile2 get their indices offset so they continue after matfile1.
          e.g. if matfile1 has 0,1 then matfile2's 0,1 -> 2,3 in outfile.
        """
        lines1 = []
        mat_indices1 = []

        # --- Read first file ---
        with open(matfile1, "r") as f1:
            for line in f1:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    lines1.append(line)  # keep comments / blank lines as-is
                    continue
                parts = stripped.split()
                try:
                    mat_id = float(parts[0])
                except ValueError:
                    # not a data line, just keep it as-is
                    lines1.append(line)
                    continue
                mat_indices1.append(mat_id)
                lines1.append(line)

        # How many materials in file 1? (assumes indices are 0,1,2,...)
        if not mat_indices1:
            raise ValueError(f"No material entries found in {matfile1}")
        offset = int(max(mat_indices1)) + 1

        # --- Process second file, remapping indices ---
        lines2 = []
        with open(matfile2, "r") as f2:
            for line in f2:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    # Optionally skip header comments from second file,
                    # or keep them with a prefix:
                    # lines2.append("# from second material file: " + stripped + "\n")
                    continue
                parts = stripped.split()
                try:
                    mat_id = float(parts[0])
                except ValueError:
                    continue  # skip weird lines

                new_id = mat_id + offset
                parts[0] = f"{new_id:.4e}"  # format like 0.0000e+00
                lines2.append(" ".join(parts) + "\n")

        # --- Write combined file ---
        with open(outfile, "w") as out:
            for line in lines1:
                out.write(line)
            for line in lines2:
                out.write(line)

        print(f"Combined material file written to {outfile}")

    def _read_material_file(self,path):
        """
        Read a material file and return:
          - comments: list of comment/blank lines (kept only from first file)
          - rows: list of (orig_id, values) where:
              orig_id is a float (material index)
              values is a list of strings for the remaining columns
        """
        comments = []
        rows = []
        with open(path, "r") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    comments.append(line)
                    continue
                parts = stripped.split()
                try:
                    mat_id = float(parts[0])
                except ValueError:
                    # Not a data line, treat as comment-ish
                    comments.append(line)
                    continue
                rows.append((mat_id, parts[1:]))
        return comments, rows

    def mapping_dict_to_array(self,mapping):
        max_old = int(max(mapping.keys()))
        arr = np.empty(max_old + 1, dtype=int)
        for k, v in mapping.items():
            arr[int(k)] = v
        return arr
       

    def combine_material_files_dedup(self,matfile1, matfile2, outfile):
        """
        Combine two material files into one, deduplicating identical materials.

        If a material in one file has identical values (all columns EXCEPT the first)
        to a material in the other file, they will share the same new index.

        Returns:
            mapping1, mapping2
            where mapping1 and mapping2 are dicts:
                original_material_index -> new_material_index
            for matfile1 and matfile2 respectively.
        """
        # Read both files
        comments1, rows1 = self._read_material_file(matfile1)
        _comments2, rows2 = self._read_material_file(matfile2)

        # Dictionary: row_key (tuple of values) -> new_id
        row_to_newid = {}
        new_rows = []  # list of (new_id, values)

        mapping1 = {}  # orig_id -> new_id for file 1
        mapping2 = {}  # orig_id -> new_id for file 2

        # Helper to register rows and build mapping
        def register_rows(rows, mapping):
            for orig_id, vals in rows:
                key = tuple(vals)
                if key in row_to_newid:
                    new_id = row_to_newid[key]
                else:
                    new_id = len(row_to_newid)
                    row_to_newid[key] = new_id
                    new_rows.append((new_id, vals))
                mapping[orig_id] = new_id

        # First register rows from file 1, then from file 2
        register_rows(rows1, mapping1)
        register_rows(rows2, mapping2)

        # Write combined file
        with open(outfile, "w") as out:
            # Write comments from first file
            for line in comments1:
                out.write(line)
            # Write unique materials, sorted by new_id
            for new_id, vals in sorted(new_rows, key=lambda x: x[0]):
                id_str = f"{float(new_id):.4e}"
                out.write(id_str + " " + " ".join(vals) + "\n")

        print(f"Combined deduplicated material file written to {outfile}")
        mapping1=self.mapping_dict_to_array(mapping1)
        mapping2=self.mapping_dict_to_array(mapping2)
        return mapping1, mapping2


