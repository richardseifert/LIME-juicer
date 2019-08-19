import os
import numpy as np

############################################
### SET THESE VARIABLES BEFORE FIRST RUN ###
############################################
# Set mainpath to your main LIME directory (the one containing src, example, etc).
mainpath = '/path/to/limemaster/'
# Set templatepath to the direcoty within mainpath that contains template mode.c files.
templatepath='templates/'

lamdict = {'CN':'cn.dat','HD':'hd.dat','CO':'co.dat','NO':'no.dat','HCN':'hcn@xpol.dat',\
    'HCO+':'hco+@xpol.dat','oH2O':'oh2o@daniel.dat','pH2O':'ph2o@daniel.dat',\
    'oNH3':'o-nh3.dat','H13CO+':'h13co+@xpol.dat','HC17O+':'hc17o+@xpol.dat',\
    'HC18O+':'hc18o+@xpol.dat','HOC+':'hco+@xpol.dat',\
    'C17O':'c17o.dat','C18O':'c18o.dat','THCO':'13co.dat','N2D+':'n2d+.dat',\
    'oH2D+':'o-h2d+.dat','pH2D+':'p-h2d+.dat','HD2+':'hd2+.dat','pH2CO':'ph2co-h2.dat',\
    'oH2CO':'oh2co-h2.dat','NO':'no_short.dat','CS':'cs@xpol.dat','pH2S':'ph2s.dat',\
    'oH2S':'oh2s.dat','SO2':'so2@xpol.dat','SO':'so@lique.dat','NO':'no_short.dat',\
    'HC15N':'hc15n@xpol.dat','H13CN':'h13cn@xpol.dat','DCO+':'dco+@xpol.dat',\
    'N2H+':'n2h+@xpol.dat','C2H':'c2h_h2.dat','N2H+':'n2h+@xpol.dat'}

isodict = {'H13CO+':'60.0','HC18O+':'500.0','oH2O':'1.333','pH2O':'4.0','oNH3':'2.0',\
    'HC17O+':'2630.0','C17O':'2630.0','C18O':'500.0','THCO':'60.0','oH2CO':'1.33',\
    'pH2CO':'4.0','oH2S':'1.33','pH2S':'4.0','HC15N':'400.0','H13CN':'60.0','HD':'33333.3'}

moldict = {'H13CO+':'HCO+','C18O':'CO','HC18O+':'HCO+','HOC+':'HCO+',\
    'oH2O':'H2O','pH2O':'H2O','ONH3':'NH3','HC17O+':'HCO+','C17O':'CO',\
    'THCO':'CO','oH2D+':'H2D+','pH2D+':'H2D+','pH2CO':'H2CO',\
    'oH2CO':'H2CO','oH2S':'H2S','pH2S':'H2S','HC15N':'HCN','H13CN':'HCN','HD':'H2','C2H':'C2H'}


def make_lime_block(trans,mol,velr,nchan,template,fits_root,\
                    distance,inclination,source_vel=0.,unit=1,\
                    pxls=250,imgres=0.06):
    '''
    ARGUMENTS:
      Required:
        trans    - int or array, Transition(s) to compute radiative transfer for.
                                 Integers should correspond to zero-indexed position of
                                 that transition in the lambda file for the molecule.

        mol      - string, Molecule to simulate emission for.

        velr     - int or array, Velocity resolution for each transition to be modeled.
                                 If single value is given, it will be used for all trans.

        nchan    - int or array, Number of channels for each transition to be modeled.
                                 If single value is given, it will be used for all trans.

        fits_root   - string, Root to be used for output .fits file.

        distance    - float, Distance to the source in pc.

        inclination - float, Inclination of source in degrees.

      Optional:
        source_vel  - float, Velocity of source.
        unit        - 0,1,2,3,or 4, Unit of the output fits cube.
                         0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
        pxls        - int, Number of pixels along each axis.
        imgres      - float, Spatial resolution in arcseconds.
    '''

    # Handle array-optional inputs
    try:
        iter(trans)
    except TypeError:
        #trans is scalar, but not for long.
        trans = [trans]
    try:
        iter(velr)
    except TypeError:
        #If scalar velr is given, vectorize it with same length as trans.
        velr = [velr for n in range(len(trans))]
    if len(velr) != len(trans):
        #Break if the wrong vector number of velr values was provided.
        raise ValueError("Must provide same number of velocity resolution "\
                         +"values as transitions")
    try:
        iter(nchan)
    except TypeError:
        #If scalar nchan is given, vectorize it with same length as trans.
        nchan = [nchan for n in range(len(trans))]
    if len(nchan) != len(trans):
        #Break if the wrong vector number of nchan values was provided.
        raise ValueError("Must provide same number of nchan "\
                         +"values as transitions")
  
    # Initialize Block.
    ntrans = len(trans)
    Block = np.zeros((ntrans,10),dtype='a120')
    
    # Create img[#] block for each transition:
    for i in range(ntrans):
        blstr = ('%d'%i).strip()
        fitsnm = fits_root+'_'+str(trans[i])+'.fits'
        Block[i][0] = "    img["+blstr+"].nchan      = %d;            // Number of channels"%(nchan[i])
        Block[i][1] = "    img["+blstr+"].velres     = %f;            // Channel resolution in m/s"%(velr[i])
        Block[i][2] = "    img["+blstr+"].trans      = %d;            // zero-indexed J quantum number"%(trans[i])
        Block[i][3] = "    img["+blstr+"].pxls       = %d;            // Pixels per dimension"%(pxls)
        Block[i][4] = "    img["+blstr+"].imgres     = %.3f;          // Resolution in arc seconds"%(imgres)
        Block[i][5] = "    img["+blstr+"].theta      = %.6f;          // 0: face-on, pi/2: edge-on"%(inclination*np.pi/180)
        Block[i][6] = "    img["+blstr+"].distance   = %.2f*PC;       // source distance in pc"%(distance)
        Block[i][7] = "    img["+blstr+"].source_vel = %.2f;          // source velocity in m/s"%(source_vel)
        Block[i][8] = "    img["+blstr+"].unit       = %d;            // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau"%(unit)
        Block[i][9] = "    img["+blstr+'].filename   = "'+fitsnm+'";  // Output filename'
   
    tmpfile = 'disktemplate_'+mol+'.c'

    f = open(mainpath+templatepath+tmpfile,'w')
    fr = open(mainpath+templatepath+template,'r')
    for linef in fr:
        if linef.strip()=='BEGINLINES':
            for b in range(ntrans):
                for t in range(10):
                    f.write(Block[b,t].decode("UTF-8")+'\n')
                f.write('\n')
        else:
            f.write(linef)
    fr.close()
    f.close()

    return tmpfile

############################################################################

def make_lime_exec(mod_in,mol,trans,velr,nchan,pint,psink,template,\
                   xroot,fits_root=None,pop_root=None,\
                   stellarmass=1.0,blending=0,lte=0,massX=1.0,abunX=1.0,\
                   deldop=100,maximumrad=600.,**blockinfo):
    '''
    This is a function that produces lime executables to run radiative transfer
    on a model disk for a given set of transitions of a single molecule.

    ARGUMENTS:
      Required:
        mod_in   - string, Path to the input model, relative to limemaster/DiskModels/.

        mol      - string, Molecule to simulate emission for.

        trans    - int or array, Transition(s) to compute radiative transfer for.
                                 Integers should correspond to zero-indexed position of
                                 that transition in the lambda file for the molecule.

        velr     - int or array, Velocity resolution for each transition to be modeled.
                                 If single value is given, it will be used for all trans.

        nchan    - int or array, Number of channels for each transition to be modeled.
                                 If single value is given, it will be used for all trans.

        pint     - int, Number of grid points to be placed within the volume of the simulation.

        psink    - int, Number of grid points to be placed across the surface of the simulation.

        template - string, Path to templace LIME model.c file, relative to limemaster/templatepath/

        xroot    - string, Root to be used for output lime.x file, as well as output .fits and
                           .pop files if separate roots are not provided.
        distance -    [in blockinfo] float, Distance to the source in pc.
        inclination - [in blockinfo] float, Inclination of source in degrees.

      Optional:
        fits_root   - string, Root to be used for output .fits file. Final form of .fits filename:
                            fitsnm = fits_root+'_'+transi+'.fits'
        pop_root    - string, Root to be used for output .pop file.
        stellarmass - float, Mass of central star in units of Msun!
        blending    - 0 or 1, LIME blending flag. (1) on, (0) off.
        lte         - 0 or 1, LIME LTE flag. (1) on, (0) off.
        massX       - float, Mass multiplier.
        abunX       - float, Abundance multiplier.
        deldop      - float, LIME doppler / turbulence parameter.
        maximumrad  - float, Maximum radius to perform RT calculations out to.
        source_vel  - [in blockinfo] float, Velocity of source.
        unit        - [in blockinfo] 0,1,2,3,or 4, Unit of the output fits cube.
                               0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
        pxls        - [in blockinfo] int, Number of pixels along each axis.
        imgres      - [in blockinfo] float, Spatial resolution in arcseconds.

    RETURNS:
        None
    '''

    #Handle dynamic defaults.
    if fits_root is None:
        fits_root = xroot
    if pop_root  is None:
        pop_root  = xroot

    #Check global dictionaries for isotope frac and lambda file.
    try:
        lamfile = lamdict[mol]
    except KeyError:
        raise ValueError('Could not identify molecule! stop.'+mol)
    try:
        isotopf = isodict[mol]
    except KeyError:
        isotopf = '1.0'

    tmpfile = make_lime_block(trans,mol,velr,nchan,template,fits_root,**blockinfo)

    #Remove some potential leftovers from previous runs.
    cmd = 'cd '+mainpath+templatepath+'; rm -f tmp1.c tmp2.c test.c'
    os.system(cmd)
    cmd = 'cd '+mainpath+'; rm -rf junk.out'
    os.system(cmd)

    #Get number of points in the input model.
    f = open('DiskModels/'+mod_in,'r')
    diskmodelpoints = len(f.readlines())
    f.close()

    #Update template 
    # Go to template directory, where lime model.c templates live.
    cmd = 'cd '+mainpath+templatepath+';'
    # Replace variables in template with values for this run.
    #   SED text editor replacement function -> 'sed s/old_text/new_text/'
    cmd +='cat %s      | sed "s/MODELFILENM/%s/"    > tmp1.c;  '%(tmpfile, mod_in.replace('/','\/'))
    cmd +='cat tmp1.c  | sed "s/POINTZINT/%i/"      > tmp2.c;  '%(pint)
    cmd +='cat tmp2.c  | sed "s/POINTZSINK/%i/"     > tmp1.c;  '%(psink)
    cmd +='cat tmp1.c  | sed "s/MASTAR/%.1f/"       > tmp2.c;  '%(stellarmass)
    cmd +='cat tmp2.c  | sed "s/DATAMOL/%s/"        > tmp1.c;  '%(lamfile)
    cmd +='cat tmp1.c  | sed "s/ISOTVAL/%s/"        > tmp2.c;  '%(isotopf)
    cmd +='cat tmp2.c  | sed "s/POPROOT/%s/"        > tmp1.c;  '%(pop_root)
    cmd +='cat tmp1.c  | sed "s/MASSX/%8.4f/"       > tmp2.c;  '%(massX)
    cmd +='cat tmp2.c  | sed "s/ABUNFRAC/%8.4f/"    > tmp1.c;  '%(abunX)
    cmd +='cat tmp1.c  | sed "s/MODFILEPT/%10i/"    > tmp2.c;  '%(diskmodelpoints)
    cmd +='cat tmp2.c  | sed "s/DOPVAL/%8.2f/"      > tmp1.c;  '%(deldop)
    cmd +='cat tmp1.c  | sed "s/BLENDFLAG/%01i/"    > tmp2.c;  '%(blending)
    cmd +='cat tmp2.c  | sed "s/MAXRAD/%10.2f/"     > tmp1.c;  '%(maximumrad)
    cmd +='cat tmp1.c  | sed "s/LTEFLAG/%01i/"      > test.c;  '%(lte)
    # Go back to main LIME directory and build LIME with new test.c
    # Then rename output executable.
    cmd +='cd '+mainpath+'; make clean; make; mv -f lime.x '+xroot+'.x;'
    # Go back to template directory and archive test.c using same root as executable.
    cmd +='cd '+mainpath+templatepath+'; cp -rf test.c '+xroot+'.c;'
    os.system(cmd)
    return

      
if __name__ == '__main__':
    import argparse
    
    #Parse arguments!!
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input_model",type=str,
                        help='Path to the input model, relative to limemaster/DiskModels/')
    parser.add_argument("-m", "--molecule"   ,type=str,
                        help='Molecule to simulate emission for.')
    parser.add_argument("-t", "--template",   type=str,
                        help='Path to templace LIME model.c file, relative to limemaster/templatepath/')
    parser.add_argument("-r", "--trans",      type=str,
                        help='Transition(s) to compute radiative transfer for, comma separated.')
    parser.add_argument("-v", "--velr",       type=str,
                        help='Velocity resolution for each transition to be modeled.')
    parser.add_argument("-n", "--nchan",      type=str,
                        help='Number of channels for each transition to be modeled.')


    parser.add_argument(      "--pint",       type=int,   default=10000,
                        help='Number of grid points to be placed within the volume of the simulation.')
    parser.add_argument(      "--psink",      type=int,   default=3000,
                        help='Number of grid points to be placed across the surface of the simulation.')

    parser.add_argument("-x", "--exec_root",  type=str,   default='lime',
                        help='Root to be used for output lime.x file')
    parser.add_argument("-f", "--fits_root",  type=str,   default=None, #Default, same as exec_root
                        help='Root to be used for output .fits file.')
    parser.add_argument("-p", "--pop_root",   type=str,   default=None, #Default, same as exec_root
                        help='Root to be used for output .pop file.')

    parser.add_argument(      "--stellarmass",type=float, default=1.0, # solar masses
                        help='Mass of central star in units of Msun')
    parser.add_argument(      "--distance",   type=float,              # parsecs
                        help='Distance to the source in pc.')
    parser.add_argument(      "--inclination",type=float, default=0.0, # degrees
                        help='Inclination of source in degrees.')
    parser.add_argument(      "--source_vel", type=float, default=0.0, # m/s
                        help='Velocity of source.')

    parser.add_argument(      "--blending",   type=int,   default=0,
                        help='LIME blending flag. (1) on, (0) off.')
    parser.add_argument(      "--lte",        type=int,   default=0,
                        help='LIME LTE flag. (1) on, (0) off.')
    parser.add_argument(      "--abunX",      type=float, default=1.0,
                        help='Mass multiplier.')
    parser.add_argument(      "--massX",      type=float, default=1.0,
                        help='Abundance multiplier.')
    parser.add_argument(      "--deldop",     type=float, default=100., # m/s
                        help='LIME doppler / turbulence parameter.')
    parser.add_argument(      "--maximumrad", type=float, default=600., # AU
                        help='Maximum radius to perform RT calculations out to.')

    parser.add_argument(      "--unit",       type=int,   default=1,
                        help='Brightness unit of the output fits cube.')
    parser.add_argument(      "--pxls",       type=int,   default=250,
                        help='Number of pixels along each axis.')
    parser.add_argument(      "--imgres",     type=float, default=0.06, # arcseconds
                        help='Spatial resolution in arcseconds.')


    #Make variables out of parsed input!!
    parsed = vars(parser.parse_args())
    mod_in    = parsed['input_model']
    mol       = parsed['molecule']
    template  = parsed['template']
    tstr      = parsed['trans']
    vstr      = parsed['velr']
    nstr      = parsed['nchan']
    pint      = parsed['pint']
    psink     = parsed['psink']
    xroot     = parsed['exec_root']
    fits_root = parsed['fits_root']
    pop_root  = parsed['pop_root']
    stellarmass=parsed['stellarmass']
    distance  = parsed['distance']
    inclination=parsed['inclination']
    source_vel= parsed['source_vel']
    blending  = parsed['blending']
    lte       = parsed['lte']
    abunX     = parsed['abunX']
    massX     = parsed['massX']
    deldop    = parsed['deldop']
    maximumrad= parsed['maximumrad']
    unit      = parsed['unit']
    pxls      = parsed['pxls']
    imgres    = parsed['imgres']


    #Parse lists of given transitions, velocity resolutions, and number of channels. 
    if ',' in tstr:
        trans = [int(n) for n in filter(None,tstr.split(','))]
    else:
        trans = int(tstr)
    if ',' in vstr:
        velr  = [int(n) for n in filter(None,vstr.split(','))]
    else:
        velr  = int(vstr)
    if ',' in nstr:
        nchan = [int(n) for n in filter(None,nstr.split(','))]
    else:
        nchan = int(nstr)

    #Generate lime executable!!
    make_lime_exec(mod_in=mod_in,mol=mol,template=template,
                  trans=trans, velr=velr,nchan=nchan,pint=pint,psink=psink,
                  xroot=xroot,fits_root=fits_root,pop_root=pop_root,
                  stellarmass=stellarmass,distance=distance,
                  inclination=inclination,source_vel=source_vel,
                  blending=blending,lte=lte,abunX=abunX,massX=massX,
                  deldop=deldop,maximumrad=maximumrad,unit=unit,
                  pxls=pxls,imgres=imgres)
