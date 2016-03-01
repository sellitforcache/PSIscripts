#! /home/l_bergmann/anaconda/bin/python -W ignore
#
# ss2dist, the MCNP surface source to histogram distribution maker
# Ryan M. Bergmann, March 2015 
# ryan.bergmann@psi.ch, ryanmbergmann@gmail.com

from pyne import mcnp
import math
import pylab, numpy, sys, cPickle, progressbar, copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import gridspec
from MCNPtools.to_energy import to_energy
from MCNPtools.to_temperature import to_temperature
from MCNPtools.to_wavelength import to_wavelength

class SourceSurf(object):
    def __init__(self):
        pass

class TrackData(object):
    def __init__(self):
        pass

class SurfSrc(mcnp.SurfSrc):
    def __init__(self, filename, mode="rb"):
        super(SurfSrc, self).__init__(filename, mode)
    def __str__(self):
        return self.print_header()
    def print_header(self):
        """Returns contents of SurfSrc's header as an informative string.

        Returns
        -------
        header_string : str
            A line-by-line listing of the contents of the SurfSrc's header.

        """
        header_string = "Code: {0} (version: {1}) [{2}]\n".format(
            self.kod, self.ver, self.loddat)
        header_string += "Problem info: ({0}) {1}\n{2}\n".format(
            self.idtm, self.probid, self.aid)
        header_string += "Showing dump #{0}\n".format(self.knod)
        header_string += (
            "{0} histories, {1} tracks, {2} record size, "
            "{3} surfaces, {4} histories\n").format(
            self.np1, self.nrss, self.ncrd,
            self.njsw, self.niss)
        header_string += (
            "{0} cells, source particle: {1},"
            " macrobody facet flag: {2}\n").format(
            self.niwr, self.mipts, self.kjaq)
        for i in self.surflist:
            header_string += (
                "Surface {0}: facet {1},"
                " type {2} with {3} parameters: (").format(
                i.id, i.facet_id, i.type, i.num_params)
            if i.num_params > 1:
                for j in i.surf_params:
                    header_string += " {0}".format(j)
            else:
                header_string += " {0}".format(i.surf_params)
            header_string += ")\n"
        header_string += "Summary Table: " + str(self.summary_table)

        return header_string
    def next_track(self):
        """Reads in track records for individual particles."""
        track_data = TrackData()
        track_info = self.get_fortran_record()
        track_data.record = track_info.get_double(abs(self.ncrd))
        track_data.nps = track_data.record[0]
        track_data.bitarray = track_data.record[1]
        track_data.wgt = track_data.record[2]
        track_data.erg = track_data.record[3]
        track_data.tme = track_data.record[4]
        track_data.x   = track_data.record[5]
        track_data.y   = track_data.record[6]
        track_data.z   = track_data.record[7]
        track_data.u   = track_data.record[8]
        track_data.v   = track_data.record[9]
        track_data.cs  = track_data.record[10]
        if ((track_data.u*track_data.u+track_data.v*track_data.v)<1.0):
            track_data.w   = math.copysign(
                math.sqrt(1 - track_data.u*track_data.u -
                          track_data.v*track_data.v), track_data.bitarray)
        else:
            track_data.w = 0.0
        return track_data
    def read_header(self):
        """Read in the header block data. This block comprises 4 fortran
        records which we refer to as: header, table1, table2, summary.
        """
        # read header record
        header = self.get_fortran_record()

        # interpret header
        self.kod = header.get_string(8)[0]  # code identifier
        if 'SF_00001' not in self.kod:
            self.ver = header.get_string(5)[0]  # code version identifier
            if '2.6.0' in self.ver:
                self.loddat = header.get_string(28)[0]  # code version date
                self.idtm = header.get_string(19)[0]    # current date and time
                self.probid = header.get_string(19)[0]  # problem id string
                self.aid = header.get_string(80)[0]     # title card of initial run
                self.knod = header.get_int()[0]        # dump number
            elif '2.7.0' in self.ver:
                self.loddat = header.get_string(28)[0]  # code version date
                self.idtm = header.get_string(19)[0]    # current date and time
                self.probid = header.get_string(19)[0]  # problem id string
                self.aid = header.get_string(80)[0]     # title card of initial run
                self.knod = header.get_long()[0]        # dump number
            elif '5    ' in self.ver:
                self.loddat = header.get_string(8)[0]  # code version date
                self.idtm = header.get_string(19)[0]    # current date and time
                self.probid = header.get_string(19)[0]  # problem id string
                self.aid = header.get_string(80)[0]     # title card of initial run
                self.knod = header.get_long()[0]        # dump number
            else:
                raise NotImplementedError("MCNP5/X Version:" +
                    self.ver.rstrip() + " not supported")

            # read table 1 record; various counts and sizes
            tablelengths = self.get_fortran_record()

            if '2.7.0' in self.ver:
            	self.np1  = tablelengths.get_long()[0]       # hist used to gen.source
            	self.nrss = tablelengths.get_long()[0]      # tracks writ. to surf.src
            	self.ncrd = tablelengths.get_long()[0]      # histories to surf.src
            	self.njsw = tablelengths.get_long()[0]      # number of surfaces
            	self.niss = tablelengths.get_long()[0]      # histories to surf.src
                #print self.np1, self.nrss, self.ncrd, self.njsw, self.niss
            	self.table1extra = list()
            	while tablelengths.num_bytes > tablelengths.pos:
            	    self.table1extra += tablelengths.get_int()	
            else:
            	# interpret table lengths
            	if '2.6.0' in self.ver:
            	    self.np1 = tablelengths.get_int()[0]    # hist used to gen. src
            	    self.nrss = tablelengths.get_int()[0]   # #tracks to surf src
            	else:
            	    self.np1 = tablelengths.get_long()[0]   # hist used to gen. src
            	    self.nrss = tablelengths.get_long()[0]  # #tracks to surf src
            	self.ncrd = tablelengths.get_int()[0]  # #values in surf src record
            	                                       # 6 for a spherical source
            	                                       # 11 otherwise
            	self.njsw = tablelengths.get_int()[0]  # number of surfaces
            	self.niss = tablelengths.get_int()[0]  # #histories to surf src
            	self.table1extra = list()
            	while tablelengths.num_bytes > tablelengths.pos:
            	    self.table1extra += tablelengths.get_int()
        if 'SF_00001' in self.kod:
            header = self.get_fortran_record()
            self.ver = header.get_string(12)[0]     # code version identifier
            self.loddat = header.get_string(9)[0]   # code version date
            self.idtm = header.get_string(19)[0]    # current date and time
            self.probid = header.get_string(19)[0]  # problem id string
            self.aid = header.get_string(80)[0]     # title card of initial run
            self.knod = header.get_int()[0]         # dump number
            # read table 1 record; various counts and sizes
            tablelengths = self.get_fortran_record()
            # interpret table lengths
            self.np1 = tablelengths.get_int()[0]     # hist used to gen.source
            self.notsure0 = tablelengths.get_int()[0]  # vals in surf src rec.
            self.nrss = tablelengths.get_int()[0]    # tracks writ. to surf.src
            self.notsure1 = tablelengths.get_int()[0]  # number of surfaces
            self.ncrd = tablelengths.get_int()[0]      # histories to surf.src
            self.njsw = tablelengths.get_int()[0]      # number of surfaces
            self.niss = tablelengths.get_int()[0]      # histories to surf.src
            self.table1extra = list()
            while tablelengths.num_bytes > tablelengths.pos:
                self.table1extra += tablelengths.get_int()
        if self.np1 < 0:
            # read table 2 record; more size info
            tablelengths = self.get_fortran_record()
            self.table2record = copy.deepcopy(tablelengths)  # copy entire record as is since values aren't changed
            self.niwr  = tablelengths.get_int()[0]   # #cells in surf.src card
            self.mipts = tablelengths.get_int()[0]   # source particle type
            self.kjaq  = tablelengths.get_int()[0]   # macrobody facet flag
            self.table2extra = list()
            while tablelengths.num_bytes > tablelengths.pos:
                self.table2extra += tablelengths.get_int()
        else:
            pass
        # Since np1 can be negative, preserve the actual np1 value while
        # taking the absolute value so that np1 can be used mathematically
        self.orignp1 = self.np1
        self.np1 = abs(self.np1)
        # get info for each surface
        self.surflist = list()
        self.surfrecordlist = list()
        if '2.7.0' in self.ver:
            for j in range(self.njsw):
                # read next surface info record
                self.surfaceinfo = self.get_fortran_record()
                self.surfrecordlist.append(self.surfaceinfo)
                surfinfo = SourceSurf()
                surfinfo.id = self.surfaceinfo.get_long()            # surface ID
                if self.kjaq == 1:
                    surfinfo.facet_id = self.surfaceinfo.get_long()  # facet ID
                else:
                    surfinfo.facet_id = -1                           # dummy facet ID
                surfinfo.type = self.surfaceinfo.get_long()           # surface type
                surfinfo.num_params = self.surfaceinfo.get_long()[0]  # #surface prm
                surfinfo.surf_params = \
                    self.surfaceinfo.get_double(surfinfo.num_params)
                self.surflist.append(surfinfo)
                #print surfinfo.id, surfinfo.facet_id, surfinfo.type, surfinfo.num_params, surfinfo.surf_params
            # we read any extra records as determined by njsw+niwr...
            # no known case of their actual utility is known currently
            for j in range(self.njsw, self.njsw+self.niwr):
                self.get_fortran_record()
                warn("Extra info in header not handled: {0}".format(j),
                     RuntimeWarning)
            # read summary table record
            summary_info = self.get_fortran_record()
            self.summaryrecord=copy.deepcopy(summary_info)
            self.summary_vec = []
            for i in range(0,summary_info.num_bytes/8):
                self.summary_vec.append(summary_info.get_long()[0])
                #print self.summary_vec[i]
            self.summary_vec=numpy.array(self.summary_vec)
            summary_info.reset()
            self.summary_table = summary_info.get_long(
                (2+4*self.mipts)*(self.njsw+self.niwr)+1)
            self.summary_extra = list()
            while summary_info.num_bytes > summary_info.pos:
                self.summary_extra += summary_info.get_long()
        else:
            for j in range(self.njsw):
                # read next surface info record
                self.surfaceinfo = self.get_fortran_record()
                self.surfrecordlist.append(self.surfaceinfo)
                surfinfo = SourceSurf()
                surfinfo.id = self.surfaceinfo.get_int()            # surface ID
                if self.kjaq == 1:
                    surfinfo.facet_id = self.surfaceinfo.get_int()  # facet ID
                else:
                    surfinfo.facet_id = -1                           # dummy facet ID
                surfinfo.type = self.surfaceinfo.get_int()           # surface type
                surfinfo.num_params = self.surfaceinfo.get_int()[0]  # #surface prm
                surfinfo.surf_params = \
                    self.surfaceinfo.get_double(surfinfo.num_params)
                self.surflist.append(surfinfo)
                print surfinfo.id, surfinfo.facet_id, surfinfo.type, surfinfo.num_params, surfinfo.surf_params
            # we read any extra records as determined by njsw+niwr...
            # no known case of their actual utility is known currently
            for j in range(self.njsw, self.njsw+self.niwr):
                self.get_fortran_record()
                warn("Extra info in header not handled: {0}".format(j),
                     RuntimeWarning)
            # read summary table record
            summary_info = self.get_fortran_record()
            self.summaryrecord=copy.deepcopy(summary_info)
            self.summary_table = summary_info.get_int(
                (2+4*self.mipts)*(self.njsw+self.niwr)+1)
            self.summary_extra = list()
            while summary_info.num_bytes > summary_info.pos:
                self.summary_extra += summary_info.get_int()

    def put_header(self):
            """Write the header part of the header to the surface source file"""
            if 'SF_00001' in self.kod:
                rec = [self.kod]
                joinrec = "".join(rec)
                newrecord = _FortranRecord(joinrec, len(joinrec))
                self.put_fortran_record(newrecord)
    
                rec = [self.ver, self.loddat, self.idtm, self.probid, self.aid]
                joinrec = "".join(rec)
                newrecord = _FortranRecord(joinrec, len(joinrec))
                newrecord.put_int([self.knod])
                self.put_fortran_record(newrecord)
            elif '2.7.0' in self.ver:
                rec = [self.kod, self.ver, self.loddat, 
                        self.idtm, self.probid, self.aid]

                joinrec = "".join(rec)
                newrecord = _FortranRecord(joinrec, len(joinrec))
                newrecord.put_long([self.knod])
                self.put_fortran_record(newrecord)
            else:
                rec = [self.kod, self.ver, self.loddat,
                       self.idtm, self.probid, self.aid]
    
                joinrec = "".join(rec)
                newrecord = _FortranRecord(joinrec, len(joinrec))
                newrecord.put_int([self.knod])
                self.put_fortran_record(newrecord)
            return

    def put_table_1(self):
        """Write the table1 part of the header to the surface source file"""
        newrecord = _FortranRecord("", 0)
        if '2.7.0' in self.ver:
            newrecord.put_long([self.orignp1 ]) 
            newrecord.put_long([self.nrss    ]) 
            newrecord.put_long([self.ncrd    ]) 
            newrecord.put_long([self.njsw    ]) 
            newrecord.put_long([self.niss    ]) 
            newrecord.put_long(self.table1extra)
            self.put_fortran_record(newrecord)
        else:
            if '2.6.0' in self.ver:
                newrecord.put_int([self.np1])
                newrecord.put_int([self.nrss])
            else:
                newrecord.put_long([self.np1])
                newrecord.put_long([self.nrss])
    
            newrecord.put_int([self.ncrd])
            newrecord.put_int([self.njsw])
            newrecord.put_int([self.niss])  # MCNP needs 'int', could be 'long' ?
            newrecord.put_int(self.table1extra)
            self.put_fortran_record(newrecord)

        return

    def put_table_2(self):
        """Write the table2 part of the header to the surface source file"""
        if '2.7.0' in self.ver:
            self.put_fortran_record(self.table2record)
        else:
            newrecord = _FortranRecord("", 0)
            newrecord.put_int([self.niwr])
            newrecord.put_int([self.mipts])
            newrecord.put_int([self.kjaq])
            newrecord.put_int(self.table2extra)
            self.put_fortran_record(newrecord)
        return

    def put_surface_info(self):
        """Write the record for each surface to the surface source file"""

        if '2.7.0' in self.ver:
            for cnt, s in enumerate(self.surflist):
                self.put_fortran_record(self.surfrecordlist[cnt])
        else:
            for cnt, s in enumerate(self.surflist):
                newrecord = _FortranRecord("", 0)
                newrecord.put_int(s.id)
                if self.kjaq == 1:
                    newrecord.put_int(s.facet_id)  # don't add a 'dummy facet ID'
                # else no macrobody flag byte in the record
    
                newrecord.put_int(s.type)
                newrecord.put_int(s.num_params)
                newrecord.put_double(s.surf_params)
    
                self.put_fortran_record(newrecord)
        return

    def put_summary(self):
        """Write the summary part of the header to the surface source file"""
        if '2.7.0' in self.ver:
            newrecord = _FortranRecord("", 0)
            newrecord.put_long(self.summary_vec)
            self.put_fortran_record(newrecord)
            #self.put_fortran_record(self.summaryrecord)
        else:
            newrecord = _FortranRecord("", 0)
            newrecord.put_int(list(self.summary_table))
            newrecord.put_int(list(self.summary_extra))
            self.put_fortran_record(newrecord)
        return

    def write_header(self):
        """Write the first part of the MCNP surface source file. The header content
        comprises five parts shown below.
        """
        self.put_header()
        self.put_table_1()
        self.put_table_2()
        self.put_surface_info()
        self.put_summary()

    def write_tracklist(self):
        """Write track records for individual particles. Second part of the MCNP
        surface source file.  Tracklist is also known as a 'phase space'.
        """

        progress = progressbar.ProgressBar()
        for j in progress(range(self.nrss)):  # nrss is the size of tracklist
            newrecord = _FortranRecord("", 0)
            # 11 records comprising particle information
            newrecord.put_double(self.tracklist[j].nps)
            newrecord.put_double(self.tracklist[j].bitarray)
            newrecord.put_double(self.tracklist[j].wgt)
            newrecord.put_double(self.tracklist[j].erg)
            newrecord.put_double(self.tracklist[j].tme)
            newrecord.put_double(self.tracklist[j].x)
            newrecord.put_double(self.tracklist[j].y)
            newrecord.put_double(self.tracklist[j].z)
            newrecord.put_double(self.tracklist[j].u)
            newrecord.put_double(self.tracklist[j].v)
            newrecord.put_double(self.tracklist[j].cs)
            self.put_fortran_record(newrecord)
        return

    def update_tracklist(self, surf_src):
        """ Update tracklist from another surface source.
        This updates the surface source in-place.
        """

        # Catch for improper non-SurfSrc type
        if type(surf_src) != SurfSrc:
            raise TypeError('Surface Source is not of type SurfSrc')

        # Because 'kod' is the first header attribute
        elif not hasattr(surf_src, 'kod'):
            raise AttributeError(
                'No header attributes for surface source argument')
        elif not hasattr(self, 'kod'):
            raise AttributeError(
                'No header attributes read for surface source')

        # Because 'tracklist' forms the non-header portion
        elif not hasattr(surf_src, 'tracklist'):
            raise AttributeError(
                'No tracklist read for surface source argument')
        elif not hasattr(self, 'tracklist'):
            raise AttributeError(
                'No tracklist read for surface source')

        # No point in updating with self
        elif self == surf_src:
            raise ValueError('Tracklist cannot be updated with itself')

        self.tracklist = surf_src.tracklist
        self.nrss = surf_src.nrss

    def __del__(self):
        """Destructor. The only thing to do is close the file."""
        self.f.close()



def make_equi_str(theta_0,N):
	import numpy
	sin2 = numpy.sin(theta_0)*numpy.sin(theta_0)
	out = [0.0]
	for i in range(0,N):
		out.append(numpy.arcsin(numpy.sqrt((i+1)*sin2)))
	return numpy.array(out)

def rotate_xy(vec,theta,deg=True):
	import numpy
	if deg:
		theta_rad = theta * numpy.pi / 180.0
	else:
		theta_rad = theta
	x=vec[0]*numpy.cos(theta_rad) - vec[1]*numpy.sin(theta_rad)
	y=vec[0]*numpy.sin(theta_rad) + vec[1]*numpy.cos(theta_rad)
	return numpy.array([x,y,vec[2]])


#
#    function to fold a fine-binned xs with a coarse spectrum
#
def rebin_xs(xs=0,xs_e=0,spec_e=0):
	# needs to be pointwise, NOT binned
	assert(len(xs_e)==len(xs))
	#print "length of xs",len(xs)
	#print "length of spec bins",len(spec_e)

	# 
	spec_xs=[]
	for i in range(0,len(spec_e)-1):
		# get requested bin widths
		low  = spec_e[i]
		high = spec_e[i+1]
		# do logic on xs E grid
		logic_low  = xs_e < low
		logic_high = xs_e > high
		dex_low  = numpy.nonzero(numpy.diff(logic_low))[0]
		dex_high = numpy.nonzero(numpy.diff(logic_high))[0]
		# figure out edge cases
		if len(dex_low) == 0:
			if logic_low[0]:	# all ones, low is above last point
				dex_low = len(xs_e)-1
			else:				# all zeros, low is below first point
				dex_low = 0
		else:
			dex_low = dex_low[0]
		if len(dex_high) == 0:
			if logic_high[0]:   # all ones, high is below first point
				dex_high = 0
			else:				# all zeros, high is above last point
				dex_high = len(xs_e)-1
		else:
			dex_high = dex_high[0]
		#print dex_low,dex_high
		# average the pointwise data 
		if dex_low == dex_high:  # bin is within two xs points
			if dex_high == len(xs_e)-1:
				e_l  = xs_e[dex_high]
				e_h  = xs_e[dex_high]
				xs_l = xs[  dex_high]
				xs_h = xs[  dex_high]
				a = 0.0
			else:
				e_l  = xs_e[dex_low]
				e_h  = xs_e[dex_high+1]
				xs_l = xs[  dex_low]
				xs_h = xs[  dex_high+1]
				a = (xs_h-xs_l)/(e_h-e_l)
			b = xs_l - a*e_l
			avg = (a/2.0)*(high*high-low*low)/(high-low)+b
		else:
			avg_vals=[]
			avg_widths=[]
			#do first bin
			e_l  = xs_e[dex_low]
			e_h  = xs_e[dex_low+1]
			xs_l = xs[  dex_low]
			xs_h = xs[  dex_low+1]
			a = (xs_h-xs_l)/(e_h-e_l)
			b = xs_l - a*e_l
			avg_vals.append( (a/2.0)*(e_h*e_h-low*low)/(e_h-low)+b )
			avg_widths.append(e_h-low)
			#do middle bins
			for i in range(dex_low,dex_high-1):
				e_l  = xs_e[i]
				e_h  = xs_e[i+1]
				xs_l = xs[  i]
				xs_h = xs[  i+1]
				a = (xs_h-xs_l)/(e_h-e_l)
				b = xs_l - a*e_l
				avg_vals.append( (a/2.0)*(e_h*e_h-e_l*e_l)/(e_h-e_l)+b )
				avg_widths.append(e_h-e_l)
			#do last bin
			if dex_high == len(xs_e)-1:
				e_l  = xs_e[dex_high]
				e_h  = xs_e[dex_high]
				xs_l = xs[  dex_high]
				xs_h = xs[  dex_high]
				a=0.0
			else:
				e_l  = xs_e[dex_high]
				e_h  = xs_e[dex_high+1]
				xs_l = xs[  dex_high]
				xs_h = xs[  dex_high+1]
				a = (xs_h-xs_l)/(e_h-e_l)
			b = xs_l - a*e_l
			avg_vals.append( (a/2.0)*(high*high-e_l*e_l)/(high-e_l)+b )
			avg_widths.append(high-e_l)
			#avg by bin width and append value
			avg_widths = numpy.array(avg_widths) # normalize
			avg = numpy.average(avg_vals,weights=avg_widths)
		spec_xs.append(avg)
	# return array
	return numpy.array(spec_xs)

def coarsen(values,bins,bin_red=2):
	import numpy
	v_out=[]
	b_out=[]
	for i in range(0,len(values)/bin_red):
		v = 0.0
		for j in range(0,bin_red):
			v = v + values[i*bin_red+j]
		v_out.append(v)
		b_out.append(bins[i*bin_red])
	b_out.append(bins[-1])
	return numpy.array(v_out),numpy.array(b_out)


def make_steps(ax,bins_in,avg_in,values_in,options=['log'],color=None,label='',ylim=False):
	import numpy, re
	assert(len(bins_in)==len(values_in)+1)

	### make copies
	bins=bins_in[:]
	values=values_in[:]
	avg=avg_in[:]
	#err=err_in[:]

	### smooth data?  parse format
	for opt in options:
		res = re.match('smooth',opt)
		if res:
			smooth_opts = opt.split('=')
			if len(smooth_opts)==1:
				wlen = 7
			elif len(smooth_opts)==2:
				wlen = int(smooth_opts[1])
			else:
				wlen = int(smooth_opts[1])
				print "MULTIPLE = SIGNS IN SMOOTH.  WHY?  ACCEPTING FIRST VALUE."
			if wlen%2==0:
				print "WINDOW LENGTH EVEN, ADDING 1..."
				wlen = wlen + 1
			print "smoothing %d bins..."%wlen
			label = label + ' SMOOTHED %d BINS'%wlen
			values = self._smooth(numpy.array(values),window_len=wlen)
			values = values[(wlen-1)/2:-(wlen-1)/2]   # trim to original length

	### coarsen data?  parse format
	for opt in options:
		res = re.match('coarsen',opt)
		print "data len ",len(values)
		if res:
			coarsen_opts = opt.split('=')
			if len(coarsen_opts)==1:
				bin_red = 2
			elif len(coarsen_opts)==2:
				bin_red = int(coarsen_opts[1])
			else:
				bin_red = int(coarsen_opts[1])
				print "MULTIPLE = SIGNS IN SMOOTH.  WHY?  ACCEPTING FIRST VALUE."
			if len(values)%bin_red==0:
				print "Reducing bins by factor of %d ..."%bin_red
				#label = label + ' COMBINED %d BINS'%bin_red
				values,bins = coarsen(numpy.array(values),numpy.array(bins),bin_red=bin_red)
			else:
				print "DATA LENGHTH NOT EVENLY DIVISIBLE BY COARSEN FACTOR, IGNORING..."


class histogram:

	def __init__(self,bins):
		self.bins		= 	copy.deepcopy(bins)  # bins are edges
		self.n_bins		=	len(bins)
		self.values		=	numpy.zeros((self.n_bins-1,))
		self.counts		= 	numpy.zeros((self.n_bins-1,))

	def add(self,bin_val,weight):

		# check if in bounds
		valid = True
		if bin_val < self.bins[0] or bin_val > self.bins[-1]:
			valid = False

		# add weight to bin if between bins
		if valid:
			dex = next((i for i, x in enumerate(bin_val < self.bins) if x), False)
			self.values[dex] = self.values[dex] + weight
			self.counts[dex] = self.counts[dex] + 1

#
#
#
#
#
#  START OF TASKS
#
#
#
#
#
charge_per_amp = 6.241e18
charge_per_milliamp = charge_per_amp/1000.0

filename = sys.argv[1]
phi_bin = int(sys.argv[2])
theta_bin = int(sys.argv[3])
#E_bin = int(sys.argv[4])
obj_bin = 1

printflag = True
errorflag = False

if printflag:
	progress = progressbar.ProgressBar()
else:
	def progress(inp):
		return inp

if filename[-4:]=='wssa':
	typeflag = 1
elif filename[-4:]=='trks':
	typeflag = 2
else:
	typeflag = 0

### interpret file name
if printflag:
	print "\n============================\n"
	if typeflag == 1:
		print "Loading '"+filename+"' as a MCNP binary surface source"
	elif typeflag == 2:
		print "Loading '"+filename+"' as a pickled pyne.mcnp.SurfSrc object"
	else:
		print "Loading '"+filename+"' as a pickled numpy.array object"

### load data
if typeflag == 1:
	ss = SurfSrc(filename)
	ss.read_header()
	#ss.read_tracklist()
	print ss.print_header()
elif typeflag == 2:
	f=open(filename,'r')
	ss = cPickle.load(f)
		
else:
	f=open(filename,'r')
	d = cPickle.load(f)
	f.close()

if printflag:
	print "Done."

this_sc = int(sys.argv[4])
type_in = sys.argv[5]
print type_in
if 'flux' in type_in:
    fluxflag=True
else:
    fluxflag=False

print fluxflag

### init
if typeflag:
	if this_sc == 81233:
		#  bin parameters
		E_bins   = numpy.array([1e-12, 1e-6, 600])
		#x_bins   = numpy.linspace(-10,10,41)
		#y_bins   = numpy.linspace(-15,25,81)
		x_bins   = numpy.array([-1.5,1.5])
		y_bins   = numpy.array([-9.0,-6.0,-3.0,0.0,3.0,6.0,9.0])
		theta_bins = numpy.array([0,90,180])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist	 = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane	= numpy.array([4.3837115E-01 , -8.9879405E-01 ,  0.0000000E+00  , 8.1346130E+02  ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([ 379.190607512  , -720.115 ,     -18.5            ])   # global again
		surface_normal 	= numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]])  
		surface_vec1    = numpy.array([-surface_plane[1], surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 81202:
		#  bin parameters
		E_bins   = numpy.array([1e-12, 1e-6, 600])
		#x_bins   = numpy.linspace(-4.75,4.75,21)
		x_bins   = numpy.linspace(-4.75,4.75,2)
		#diff     = x_bins[1]-x_bins[0]
		#x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		#x_bins   = numpy.append(x_bins,  x_bins[-1]+diff)
		#y_bins   = numpy.linspace(-8.5,8.5,41)
		y_bins   = numpy.linspace(-8.5,8.5,2)
		#diff     = y_bins[1]-y_bins[0]
		#y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		#y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		theta_bins = numpy.array([0,90,180])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist	 = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters 
		surface_plane	= numpy.array([4.3837115E-01 , -8.9879405E-01 ,  0.0000000E+00 ,  2.1106130E+02  ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([115.105, -178.686751185 ,     -18.5            ])   # global again
		surface_normal 	= numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]])  
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])   # why negative?! - because normal is OUTWARDS away from target
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
	elif this_sc == 20359:
		#  bin parameters
		E_bins   = numpy.array([1e-12, 1e-6 ,600])
		#x_bins   = numpy.linspace(-7,7,41)
		x_bins   = numpy.linspace(-7,7,29)
		#diff     = x_bins[1]-x_bins[0]
		#x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		#x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		#y_bins   = numpy.linspace(-13.5,13.5,81)
		y_bins   = numpy.linspace(-13.5,13.5,55)
		#diff     = y_bins[1]-y_bins[0]
		#y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		#y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		#E_bins	 = to_energy(numpy.array([0,0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,10]))
		#E_bins   = E_bins[::-1]
		#x_bins   = numpy.linspace(-7,7,5)
		#y_bins   = numpy.linspace(-13.5,13.5,5)
		theta_bins = numpy.array([0,90,180])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist	 = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		#surface_plane	= numpy.array([-4.0673664E-01 ,  9.1354546E-01 ,  0.0000000E+00 , -3.0003472E+00  ])   # plane, GLOBAL coordinates
		surface_plane	= numpy.array([-4.0673664E-01 ,  9.1354546E-01 ,  0.0000000E+00 , -5.0003472E+00])
		surface_center  = numpy.array([ 24.095 ,  7.443876237  ,     -16.            ])   # global again
		surface_normal 	= numpy.array([-surface_plane[0],-surface_plane[1],-surface_plane[2]])   # SS written for -20359!  normal is flipped
		surface_vec1    = numpy.array([ surface_plane[1],-surface_plane[0] ,  0.0])               # why negative?! - because normal is OUTWARDS away from target
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
	elif this_sc == 10146:
		#  bin parameters
		E_bins   = numpy.array([1e-12, 1e-6 ,60])
		#x_bins   = numpy.linspace(-7,7,21)
		x_bins   = numpy.linspace(-2.5,2.5,11)
		diff     = x_bins[1]-x_bins[0]
		x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		#y_bins   = numpy.linspace(-13.5,13.5,41)
		y_bins   = numpy.linspace(-6,6,25)
		diff     = y_bins[1]-y_bins[0]
		y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		#E_bins  = to_energy(numpy.array([0,0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,10]))
		#E_bins   = E_bins[::-1]
		#x_bins   = numpy.linspace(-7,7,5)
		#y_bins   = numpy.linspace(-13.5,13.5,5)
		theta_bins = numpy.array([0,90,180])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane   = numpy.array([0.99863,  -0.052336 , 0.0 , 623.0  ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([ 623.5654 ,  -5.52  ,    0.            ])   # global again
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 4:
		#  bin parameters
		E_bins   = numpy.array([1e-12, 1e-6 ,60])
		x_bins   = numpy.linspace(-15,15,81)
		diff     = x_bins[1]-x_bins[0]
		x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		y_bins   = numpy.linspace(-15,15,81)
		diff     = y_bins[1]-y_bins[0]
		y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		#E_bins  = to_energy(numpy.array([0,0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,10]))
		#E_bins   = E_bins[::-1]
		#x_bins   = numpy.linspace(-7,7,5)
		#y_bins   = numpy.linspace(-13.5,13.5,5)
		theta_bins = numpy.array([0,90,180])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane   = numpy.array([   1,  0.0 , 0.0 , 148  ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([ 148 ,  0  ,    0.            ])   # global again
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_vec1    = numpy.array([surface_plane[1], -surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 10150:
		#  bin parameters
		E_bins   = numpy.array([0,  to_energy(3)])  # 2.27e-9 = 6 A, 9.09e-9 = 3 A ,1e-6 = 0.3 A
		#E_bins   = to_energy(numpy.linspace(1,13,121))
		#E_bins   = E_bins[::-1]
		#x_bins   = numpy.linspace(-1,2,4)
		x_bins   = numpy.linspace(-15,15,61)
		diff     = x_bins[1]-x_bins[0]
		#x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		#x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		#y_bins   = numpy.linspace(-4,1,6)
		y_bins   = numpy.linspace(-15,15,61)
		diff     = y_bins[1]-y_bins[0]
		#y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		#y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		#E_bins  = to_energy(numpy.array([0,0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,10]))
		#E_bins   = E_bins[::-1]
		#x_bins   = numpy.linspace(-7,7,5)
		#y_bins   = numpy.linspace(-13.5,13.5,5)
		#theta_bins = numpy.array([0,0.25,0.5,0.75,1.0,1.25,1.50,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,90,180])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		theta_bins = numpy.array([0,2.0])*numpy.pi/180.0
		#theta_bins  = make_equi_str(0.5*numpy.pi/180.0,64)
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		# surface_plane   = numpy.array([   0.99961, 0.0279216, 0.0, 159.543  ])   # plane, GLOBAL coordinates
		surface_plane   = numpy.array([  0.99863 ,0.052336 ,0.0 , 159.543  ])
		#surface_center  = numpy.array([ 160.5432 ,  -33.58  ,    0.            ])   # global again
		# centered on AMOR guide
		surface_center  = numpy.array([ 160.436237733 ,  -29.75  ,    0.            ])
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 10156:
		#  bin parameters
		E_bins   = numpy.array([1e-12, 1e-6 ,600])
		x_bins   = numpy.linspace(-25,25,81)
		diff     = x_bins[1]-x_bins[0]
		x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		y_bins   = numpy.linspace(-20,20,81)
		diff     = y_bins[1]-y_bins[0]
		y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		#E_bins  = to_energy(numpy.array([0,0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,10]))
		#E_bins   = E_bins[::-1]
		#x_bins   = numpy.linspace(-7,7,5)
		#y_bins   = numpy.linspace(-13.5,13.5,5)
		theta_bins = numpy.array([0,90,180])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane   = numpy.array([  0.990268 , -0.1391730 , 0.0, 165.728  ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([ 159.0746496,  -58.93  ,    0.            ])   # global again
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 10110:
		#  bin parameters
		E_bins   = numpy.array([1e-12, 5e-9, 1e-6 ,600])  # 5e-9 = 4 A, 1e-6 = 0.3 A
		x_bins   = numpy.linspace(-2.5,2.5,11)
		diff     = x_bins[1]-x_bins[0]
		x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		y_bins   = numpy.linspace(-6,6,25)
		diff     = y_bins[1]-y_bins[0]
		y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		#E_bins  = to_energy(numpy.array([0,0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,10]))
		#E_bins   = E_bins[::-1]
		#x_bins   = numpy.linspace(-7,7,5)
		#y_bins   = numpy.linspace(-13.5,13.5,5)
		theta_bins = numpy.array([0,0.125,2,90,180])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane   = numpy.array([   0.99863 ,  0.052336 , 0.0 , 159.573 ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([  161.3210058 , -29.75  ,    0.            ])   # global again
		#surface_plane   = numpy.array([   0.99961, 0.0279216, 0.0, 159.573 ])   # plane, GLOBAL coordinates
		#surface_center  = numpy.array([ 160.5432 ,  -33.58  ,    0.            ])   # global again
		# centered on AMOR guide
		#surface_center  = numpy.array([ 160.436237733 ,  -29.75  ,    0.            ])
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 10160:
		#  bin parameters
		E_bins   = numpy.array([1e-12,600])  # 5e-9 = 4 A, 1e-6 = 0.3 A
		x_bins   = numpy.linspace(-2.0,2.0,11)
		diff     = x_bins[1]-x_bins[0]
		x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		y_bins   = numpy.linspace(-6,6,25)
		diff     = y_bins[1]-y_bins[0]
		y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		#E_bins  = to_energy(numpy.array([0,0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,10]))
		#E_bins   = E_bins[::-1]
		#x_bins   = numpy.linspace(-7,7,5)
		#y_bins   = numpy.linspace(-13.5,13.5,5)
		theta_bins = numpy.array([0,0.125,2,90,180])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane   = numpy.array([   0.999994 , 0.00349065 , 0.0 , 161.0 ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([  161.132250181, -37.61  ,    0.            ])   # global again
		#surface_plane   = numpy.array([   0.99961, 0.0279216, 0.0, 159.573 ])   # plane, GLOBAL coordinates
		#surface_center  = numpy.array([ 160.5432 ,  -33.58  ,    0.            ])   # global again
		# centered on AMOR guide
		#surface_center  = numpy.array([ 160.436237733 ,  -29.75  ,    0.            ])
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 10010:
		#  bin parameters
		E_bins   = numpy.array([1e-12,9.09e-9,1e-6,600])  # 5e-9 = 4 A, 9.09e-9 = 3A, 1e-6 = 0.3 A
		x_bins   = numpy.linspace(-12,12,51)
		#diff     = x_bins[1]-x_bins[0]
		#x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		#x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		y_bins   = numpy.linspace(-12,12,51)
		#diff     = y_bins[1]-y_bins[0]
		#y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		#y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		theta_bins = numpy.array([0,90])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		#theta_bins  = make_equi_str(10.0*numpy.pi/180.0,32)
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane   = numpy.array([   9.9862953E-01 , -5.2335956E-02  , 0.0000000E+00  , 1.6486150E+01 ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([  14.48, -38.71  ,    0.            ])   # global again
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_normal_rot = surface_normal 
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 10171:
		#  bin parameters
		E_bins   = numpy.array([1e-12,9.09e-9,1e-6,600])  # 5e-9 = 4 A, 9.09e-9 = 3A, 1e-6 = 0.3 A
		#E_bins    = numpy.array([8.5124057517656764e-09,9.7270177496394964e-09])
		x_bins   = numpy.linspace(-12,12,49)
		#diff     = x_bins[1]-x_bins[0]
		#x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		#x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		y_bins   = numpy.linspace(-12,12,49)
		#diff     = y_bins[1]-y_bins[0]
		#y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		#y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		#theta_bins = numpy.array([0,0.25,0.5,0.75,1.0,180])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		theta_bins = numpy.array([0,90])*numpy.pi/180.0 
		#theta_bins  = make_equi_str(1.0*numpy.pi/180.0,10)
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane   = numpy.array([  0.99863 , -0.052336 , 0.0 , 130.2001 ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([  128.5874203,  -34.18,   0.  ])   # global again
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_normal_rot = rotate_xy(surface_normal,3.0)
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 10172:
		#  bin parameters
		E_bins   = numpy.array([1e-12,1e-6,600])
		#x_bins   = numpy.linspace(-45,45,181)
		#x_bins   = numpy.linspace(-20,20,41)
		#x_bins   = numpy.linspace(-6.76,-1.76,11)
		x_bins   = numpy.linspace(-11,11,45)
		#diff     = x_bins[1]-x_bins[0]
		#x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		#x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		#y_bins   = numpy.linspace(-45,45,181)
		#y_bins   = numpy.linspace(-10,10,21)
		y_bins   = numpy.linspace(-7,7,29)
		#theta_bins = numpy.array([0,0.25,0.5,0.75,1.0,180])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		#theta_bins = numpy.array([0.00,0.15,0.30,0.60,1.00,1.50,2.0,2.5,3.0])*numpy.pi/180.0 
		theta_bins = numpy.array([0.00,0.60,1.00,1.50,2.0,2.5,3.0])*numpy.pi/180.0 
		#theta_bins  = make_equi_str(1.0*numpy.pi/180.0,10)
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane   = numpy.array([ 0.99863 , -0.052336 , 0.0 , 130.2002])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([  127.5438191,  -54.095,   0.  ])   # global again
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		yz_rotation_degrees =  0.0
		xy_rotation_degrees = -6.0
		surface_normal_rot  = rotate_xy(surface_normal,xy_rotation_degrees)  # FOCUS is -6 degrees off plane
		surface_vec1_rot    = rotate_xy(surface_vec1,  xy_rotation_degrees)  # FOCUS is -6 degrees off plane
		surface_vec2_rot    = rotate_xy(surface_vec2,  xy_rotation_degrees)  # FOCUS is -6 degrees off plane
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 10177:
		#  bin parameters
		E_bins   = numpy.array([1e-12,9.09e-9,1e-6,600])  # 5e-9 = 4 A, 9.09e-9 = 3A, 1e-6 = 0.3 A
		x_bins   = numpy.linspace(-3,3,13)
		diff     = x_bins[1]-x_bins[0]
		x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		y_bins   = numpy.linspace(-7,7,29)
		diff     = y_bins[1]-y_bins[0]
		y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		#theta_bins = numpy.array([0,0.25,0.5,0.75,1.0,180])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		#theta_bins = numpy.array([0,90])*numpy.pi/180.0 
		theta_bins  = make_equi_str(1.0*numpy.pi/180.0,3)
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane   = numpy.array([  0.999994  , 0.00349065  ,  0.0 , 160.604 ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([  160.7362478 , -37.61 ,   0.  ])   # global again
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_normal_rot = surface_normal
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 10113:
		#  bin parameters
		#expon = numpy.linspace(-11,3,1025)
		#E_bins   =   numpy.power(10.0,expon) 
		#E_bins   = numpy.array([1e-12,1e-6,10,600])  # 5e-9 = 4 A, 9.09e-9 = 3A, 1e-6 = 0.3 A
		E_bins   = numpy.array([1e-6,600])
		#x_bins   = numpy.linspace(-45,45,181)
		#x_bins   = numpy.linspace(-20,20,41)
		x_bins   = numpy.linspace(-19.1,-14.1,11)
		#diff     = x_bins[1]-x_bins[0]
		#x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		#x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		#y_bins   = numpy.linspace(-45,45,181)
		#y_bins   = numpy.linspace(-10,10,21)
		y_bins   = numpy.linspace(-7,7,29)
		#diff     = y_bins[1]-y_bins[0]
		#y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		#y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		theta_bins = numpy.array([0,5.0])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		#theta_bins  = make_equi_str(10.0*numpy.pi/180.0,32)
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane   = numpy.array([  0.99863,  -0.052336,  0.0,  130.11 ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([  128.007636,  -43.52137 ,    0. ])
		#surface_center  = numpy.array([  128.059972,  -42.52274 ,    0.      ])
		#surface_center  = numpy.array([  127.9553,  -44.52,  0  ])   # global again
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_normal_rot = surface_normal 
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=1024.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 10175:
		#  bin parameters
		#expon = numpy.linspace(-11,3,1025)
		#E_bins   =   numpy.power(10.0,expon) 
		E_bins   = numpy.array([1e-12,2.27e-9,9.09e-9,1e-6,600])  # 5e-9 = 4 A, 9.09e-9 = 3A, 1e-6 = 0.3 A, 2.27e-9 = 6 A, 9.09e-9 = 3 A ,1e-6 = 0.3 A
		x_bins   = numpy.linspace(-3,3,25)
		#diff     = x_bins[1]-x_bins[0]
		#x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		#x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		y_bins   = numpy.linspace(-7,7,57)
		#diff     = y_bins[1]-y_bins[0]
		#y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		#y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		theta_bins = numpy.array([0,90])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		#theta_bins  = make_equi_str(10.0*numpy.pi/180.0,32)
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane   = numpy.array([  0.99863,  0.052336,  0.0,  158.746 ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([  160.5229124, -29.75, 0.0 ])
		#surface_center  = numpy.array([  128.059972,  -42.52274 ,    0.      ])
		#surface_center  = numpy.array([  127.9553,  -44.52,  0  ])   # global again
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_normal_rot = surface_normal 
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
	elif this_sc == 10189:
		#  bin parameters
		#expon = numpy.linspace(-11,3,1025)
		#E_bins   =   numpy.power(10.0,expon) 
		E_bins   = numpy.array([1e-12,2.27e-9,9.09e-9,1e-6,600])  # 5e-9 = 4 A, 9.09e-9 = 3A, 1e-6 = 0.3 A, 2.27e-9 = 6 A, 9.09e-9 = 3 A ,1e-6 = 0.3 A
		x_bins   = numpy.linspace(-3,3,25)
		#diff     = x_bins[1]-x_bins[0]
		#x_bins   = numpy.insert(x_bins,0,x_bins[0] -diff)
		#x_bins   = numpy.append(x_bins,  x_bins[-1]+diff) 
		y_bins   = numpy.linspace(-7,7,57)
		#diff     = y_bins[1]-y_bins[0]
		#y_bins   = numpy.insert(y_bins,0,y_bins[0] -diff)
		#y_bins   = numpy.append(y_bins,  y_bins[-1]+diff)
		theta_bins = numpy.array([0,90])*numpy.pi/180.0   # 90 included as sanity check, ss should only write tracks in normal dir
		#theta_bins  = make_equi_str(10.0*numpy.pi/180.0,32)
		phi_bins = numpy.linspace(0,2*numpy.pi,2) 
		dist     = numpy.zeros((  len(E_bins)-1 , len(theta_bins)-1 , len(phi_bins)-1 , len(y_bins)-1 , len(x_bins)-1 ),dtype=numpy.float64)
		#  surface plane parameters
		surface_plane   = numpy.array([  0.99863,  0.052336,  0.0,  622.746 ])   # plane, GLOBAL coordinates
		surface_center  = numpy.array([  623.8867411, -5.465, 0.0 ])
		#surface_center  = numpy.array([  128.059972,  -42.52274 ,    0.      ])
		#surface_center  = numpy.array([  127.9553,  -44.52,  0  ])   # global again
		surface_normal  = numpy.array([surface_plane[0],surface_plane[1],surface_plane[2]]) 
		surface_normal_rot = surface_normal 
		surface_vec1    = numpy.array([-surface_plane[1],surface_plane[0] ,  0.0])
		surface_vec2    = numpy.array([0.0,0.0,1.0])
		# spectrum plot
		spec_res=256.
		surface_area=(x_bins[-1]-x_bins[0])*(y_bins[-1]-y_bins[0])
else:
	dist            	= d['dist']           
	E_bins          	= d['E_bins']         
	x_bins          	= d['x_bins']         
	y_bins          	= d['y_bins']         
	theta_bins      	= d['theta_bins']     
	phi_bins        	= d['phi_bins']       
	surface_plane   	= d['surface_plane']  
	surface_normal  	= d['surface_normal']
	surface_center  	= d['surface_center'] 
	#surface_center  = numpy.array([  127.5438191,  -54.095,   0.  ])   # global again
	surface_vec1    	= d['surface_vec1']   
	surface_vec2    	= d['surface_vec2']
	surface_normal_rot	= d['surface_normal_rot'] 
	surface_vec1_rot	= d['surface_vec1_rot']   
	surface_vec2_rot	= d['surface_vec2_rot']   
	xy_rotation_degrees	= d['xy_rotation_degrees']
	yz_rotation_degrees	= d['yz_rotation_degrees']
	surface_nps     	= d['surface_nps']    
	total_weight    	= d['total_weight']   
	total_tracks    	= d['total_tracks']   
	npstrack_ratio  	= d['npstrack_ratio'] 
	this_sc         	= d['this_sc']            
	histograms_curr		= d['histograms_curr']
	histograms_flux		= d['histograms_flux']
	spec_res			= d['spec_res']

### print some details
if printflag:
	print "\n============================\n"

	if typeflag == 1:
		print "Binning '"+filename+"' according to:\n"
	else:
		print "Binning in '"+filename+"' done according to:\n"

	print "***NEUTRONS ONLY***"
	print ""
	print "Energy bin boundaries (MeV)\n",E_bins 
	print ""
	print "Wavelength bin boundaries (A)\n",to_wavelength(E_bins)
	print "    "
	print "Theta (polar) bin boundaries (degrees)\n", theta_bins*180.0/numpy.pi
	print "    "
	print "Phi (azimuthal) bin boundaries (degrees)\n", phi_bins*180.0/numpy.pi
	print "    "
	print "Y bin boundaries (cm)\n", y_bins
	print "    "
	print "X bin boundaries (cm)\n", x_bins
	print "    "
	print "NORMAL HAS BEEN ROTATED (only for angular calculations, not position!):\n",
	print "   X-Y:  % 4.2f degrees" % xy_rotation_degrees
	print "   Y-Z:  % 4.2f degrees" % yz_rotation_degrees
	print "    "

### check to make sure surface basis vectors are orthogonal
assert( numpy.abs(numpy.dot(surface_vec1,surface_vec2)  ) <= 1.e-8 )
assert( numpy.abs(numpy.dot(surface_vec1,surface_normal)) <= 1.e-8 )
assert( numpy.abs(numpy.dot(surface_vec2,surface_normal)) <= 1.e-8 )

### plot positions and vectors to make sure everything is OK

# 3d plot objects
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# origin
ax.scatter([0.0],[0.0],[0.0],'o',color='r')
# center of plane
ax.scatter(surface_center[0],surface_center[1],surface_center[2],'o',color='b')
# plane lines
plane_size = 0.3*numpy.linalg.norm(surface_center)
p1 = surface_center + surface_vec2*plane_size - surface_vec1*plane_size
p2 = surface_center + surface_vec2*plane_size + surface_vec1*plane_size
p3 = surface_center - surface_vec2*plane_size + surface_vec1*plane_size
p4 = surface_center - surface_vec2*plane_size - surface_vec1*plane_size
pts=numpy.vstack((p1,p2,p3,p4,p1))
ax.plot(pts[:,0],pts[:,1],pts[:,2],color='b')
# 
# print surface_normal
# print surface_plane[0:3]
# print numpy.cross(surface_plane[0:3],surface_normal)
# print pts
#
# normal
ax.quiver(surface_center[0],surface_center[1],surface_center[2],surface_normal[0],surface_normal[1],surface_normal[2],            color='b',pivot='tail',length=plane_size)
ax.quiver(surface_center[0],surface_center[1],surface_center[2],  surface_vec1[0],  surface_vec1[1],  surface_vec1[2],            color='b',pivot='tail',length=plane_size)
ax.quiver(surface_center[0],surface_center[1],surface_center[2],  surface_vec2[0],  surface_vec2[1],  surface_vec2[2],            color='b',pivot='tail',length=plane_size)
# rotated normal
ax.quiver(surface_center[0],surface_center[1],surface_center[2],surface_normal_rot[0],surface_normal_rot[1],surface_normal_rot[2],color='r',pivot='tail',length=plane_size)
ax.quiver(surface_center[0],surface_center[1],surface_center[2],  surface_vec1_rot[0],  surface_vec1_rot[1],  surface_vec1_rot[2],color='r',pivot='tail',length=plane_size)
ax.quiver(surface_center[0],surface_center[1],surface_center[2],  surface_vec2_rot[0],  surface_vec2_rot[1],  surface_vec2_rot[2],color='r',pivot='tail',length=plane_size)
# show 
plt.show()

### make bins for histograms
finebins=[]
Emin=E_bins[0]
Emax=E_bins[-1]
for j in range(0,int(spec_res)+1):
	finebins.append(Emin*numpy.power(Emax/Emin, j/spec_res))
finebins = numpy.array(finebins)
avg=(finebins[:-1]+finebins[1:])/2.0

### scan tracks
if typeflag:

	### init some stuff
	last_nps = 0
	x_avg = 0.
	x_dex_avg = 0
	bitarrays=[]
	histograms_curr=[]
	histograms_flux=[]
	wgt_avg = 0.5
	count =1.0
	for i in range(0,len(theta_bins)-1):
		histograms_curr.append(histogram(finebins))
		histograms_flux.append(histogram(finebins))

	if printflag:
		print "\n============================\n"
		print "Binning tracks... "
	for i in progress(range(1,ss.nrss)):
		
		### get track global position/direction
		track = ss.next_track()

		### decode bitarray
		b   = abs(track.bitarray)      # sign means what?
		j   = int(b / 2e8)             # collided?  history?
		ipt = int(b / 1e6 - j*2e2)     # particle type (1=n,2=p,3=e,4=mu-,9=proton,20=pi_+)
		nsf = int(b - ipt*1e6 - j*2e8) # surface
		
		### get data
		vec = numpy.array([track.u,track.v,track.w])
		pos = numpy.array([track.x,track.y,track.z])
		this_E 	  = track.erg
		this_wgt  = track.wgt
		
		### mcnp6
		if 'SF_00001' in ss.kod:
			nsf=track.cs
			ipt=1  #have manually set, IS NOT READ HERE, SCRIPT WILL ASSUME ALL ARE NEUTRONS
		sense = surface_plane[0]*pos[0] + surface_plane[1]*pos[1] + surface_plane[2]*pos[2] - surface_plane[3]  # use sense almost zero for on-plane particles since I don't think mcnpx prints which surface its on!

		if  (ipt==1) and (abs(sense)<=1e-5): # (nsf==this_sc): #

			### transform vector to normal system
			this_vec = numpy.array([numpy.dot(surface_vec1,vec),numpy.dot(surface_vec2,vec),numpy.dot(surface_normal_rot,vec)])
			#this_vec = numpy.array([numpy.dot(surface_vec1,vec),numpy.dot(surface_vec2,vec),numpy.dot(surface_normal,vec)])

			### transform position to surface coordinates using basis vectors specified
			xfm_pos  = numpy.subtract(pos,surface_center)
			this_pos = numpy.array([numpy.dot(surface_vec1,xfm_pos),numpy.dot(surface_vec2,xfm_pos)])
		
			### calc angular values
			this_theta  = numpy.arccos(this_vec[2])
			this_phi = numpy.arctan2(this_vec[1],this_vec[0])
			if this_phi < 0.0:
				this_phi = 2.0*numpy.pi + this_phi
		
			### find the bins
			if (this_E > E_bins[0] and this_E < E_bins[-1]):
				E_dex 	=  numpy.nonzero(this_E      < E_bins  )[0][0]-1
			else:
				E_dex = sys.maxint
			if (this_pos[0] > x_bins[0] and this_pos[0] < x_bins[-1]):
				x_dex 	=  numpy.nonzero(this_pos[0] < x_bins  )[0][0]-1
			else:
				x_dex= sys.maxint
			if (this_pos[1] > y_bins[0] and this_pos[1] < y_bins[-1]):	
				y_dex 	=  numpy.nonzero(this_pos[1] < y_bins  )[0][0]-1
			else:
				y_dex= sys.maxint
			if (this_theta > theta_bins[0] and this_theta < theta_bins[-1]):
				theta_dex	=  numpy.nonzero(this_theta     < theta_bins )[0][0]-1
			else:
				theta_dex= sys.maxint
			if (this_phi > phi_bins[0] and this_phi < phi_bins[-1]):	
				phi_dex	=  numpy.nonzero(this_phi    < phi_bins)[0][0]-1
			else:
				phi_dex=sys.maxint
				
			### increment array
			if (E_dex < len(E_bins)-1) and (theta_dex < len(theta_bins)-1) and (phi_dex < len(phi_bins)-1) and (y_dex < len(y_bins)-1) and (x_dex < len(x_bins)-1) :
				count = count+1
				wgt_avg = wgt_avg*(count-1.0)/count + this_wgt/count
				if this_wgt <= 20.0*wgt_avg:
					#if this_wgt>wgt_avg:
					#	wgt_avg=this_wgt
					x_avg = x_avg + x_bins[x_dex]
					x_dex_avg = x_dex_avg + x_dex
					#if this_wgt<=1.0:
					if fluxflag:
						dist[E_dex][theta_dex][phi_dex][y_dex][x_dex] = dist[E_dex][theta_dex][phi_dex][y_dex][x_dex] + this_wgt/this_vec[2]
						cosines[theta_dex].append(this_vec[2])
					else:
						dist[E_dex][theta_dex][phi_dex][y_dex][x_dex] = dist[E_dex][theta_dex][phi_dex][y_dex][x_dex] + this_wgt
					histograms_curr[theta_dex].add(this_E,this_wgt)
					histograms_flux[theta_dex].add(this_E,this_wgt/this_vec[2]/surface_area)
				#if this_wgt < 0:
				#	print this_wgt
				#print "accepted",this_E,x_dex,y_dex,this_wgt
				#print this_E," between ", E_bins[E_dex:E_dex+2]
			else:
				if (E_dex >= len(E_bins)-1 and printflag and errorflag): 
					print "E = %6.4E index %i is outside bin boundaries" % (this_E,E_dex,)
				if(theta_dex >= len(theta_bins)-1 and printflag and errorflag): 
					print "theta = %6.4E index %i is outside bin boundaries" % (this_theta,theta_dex)
				if(phi_dex >= len(phi_bins)-1 and printflag and errorflag): 
					print "phi = %6.4E index %i is outside bin boundaries" % (this_phi,phi_dex)
				if(y_dex >= len(y_bins)-1 and printflag and errorflag): 
					print "y = %6.4E index %i is outside bin boundaries" % (this_pos[1],y_dex)
				if(x_dex >= len(x_bins)-1 and printflag and errorflag):
					print "x = %6.4E index %i is outside bin boundaries" % (this_pos[0],x_dex)
	print "max weight",wgt_avg
	### normalize dist to nps:
	unit_area = (y_bins[1]-y_bins[0])*(x_bins[1]-x_bins[0])
	surface_nps = abs(track.nps)
	total_weight = 0.0
	total_tracks = 0
	for i in range(0,len(theta_bins)-1):
		total_tracks = total_tracks + numpy.sum(histograms_curr[i].counts)
		total_weight = total_weight + numpy.sum(histograms_curr[i].values)
		histograms_curr[i].values = histograms_curr[i].values / surface_nps
		histograms_flux[i].values = histograms_flux[i].values / surface_nps
	npstrack_ratio = surface_nps/total_tracks
	if fluxflag:
		dist = dist / surface_nps / unit_area
	else:
		dist = dist / surface_nps


	### dump array to file
	if printflag:
		print "\n============================\n"
		print "writing binned array to file 'dist'... "
	f = open('dist','wf')
	d = {}
	d['dist']=dist
	d['E_bins']=E_bins
	d['x_bins']=x_bins
	d['y_bins']=y_bins
	d['theta_bins']=theta_bins
	d['phi_bins']=phi_bins
	d['surface_plane']=surface_plane
	d['surface_center']=surface_center
	d['surface_normal']=surface_normal
	d['surface_vec1']=surface_vec1
	d['surface_vec2']=surface_vec2
	d['surface_normal_rot']=surface_normal_rot
	d['surface_vec1_rot']=surface_vec1_rot
	d['surface_vec2_rot']=surface_vec2_rot
	d['xy_rotation_degrees']=xy_rotation_degrees
	d['yz_rotation_degrees']=yz_rotation_degrees
	d['surface_nps']=surface_nps
	d['total_weight']=total_weight
	d['total_tracks']=total_tracks
	d['npstrack_ratio']=npstrack_ratio
	d['this_sc']=this_sc
	d['histograms_curr']=histograms_curr
	d['histograms_flux']=histograms_flux
	d['spec_res']=spec_res

	cPickle.dump(d,f)
	f.flush()
	f.close()
	if printflag:
		print "Done."

if printflag:
	print "\n============================\n"

if typeflag == 1:
	ss.close()


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=16)


### images
zap_x1=[-6.6, -19.1, -19.1, -6.6, -6.6]
zap_x2=[4.6,19.1 ,19.1 ,4.6, 4.6]
zap_y=[7.05, 7.05 ,-7.05, -7.05, 7.05]
x_AMOR=[2.5,2.5,-2.5,-2.5,2.5]
y_AMOR=[-6,6,6,-6,-6]
x_FOCUS=[-1.76,-1.76,-6.76,-6.76,-1.76]
y_FOCUS=[-6,6,6,-6,-6]
for theta_bin in range(0,len(theta_bins)-1):
	for E_bin in range(0,len(E_bins)-1):
		f = plt.figure()
		ax = f.add_subplot(111)
		imgplot = ax.imshow(dist[E_bin][theta_bin][phi_bin][:][:],extent=[x_bins[0],x_bins[-1],y_bins[0],y_bins[-1]],origin='lower',cmap=plt.get_cmap('jet'))
		this_weight = numpy.sum(dist[E_bin][theta_bin][phi_bin][:][:])#/surface_nps
		imgplot.set_interpolation('nearest')
		theta_deg = theta_bins[theta_bin:theta_bin+2]*180.0/numpy.pi
		phi_deg = phi_bins[phi_bin:phi_bin+2]*180.0/numpy.pi
		E_meV   = E_bins[E_bin:E_bin+2]*1.0e9
		E_eV   = E_bins[E_bin:E_bin+2]*1.0e6
		ax.set_ylabel(r'y (cm)')
		ax.set_xlabel(r'x (cm)')
		ax.plot(x_FOCUS,y_FOCUS,'0.5',linewidth=4,linestyle='--')
		#ax.set_title(r'Energies %4.2f - %4.2f meV \\       $\theta$ %4.2f - %4.2f $^{\circ}$, $\phi$ %4.2f - %4.2f $^{\circ}$ \\ nps %d tracks %d \\ total weight/nps %4.2E' % (E_meV[0],E_meV[1],theta_deg[0],theta_deg[1],phi_deg[0],phi_deg[1],int(surface_nps),int(track_count[E_bin]),this_weight))
		ax.set_title(r'Energies %4.2E - %4.2E eV \\       $\theta$ %4.2f - %4.2f $^{\circ}$, $\phi$ %4.2f - %4.2f $^{\circ}$ \\ Total weight/nps %4.2E' % (E_eV[0],E_eV[1],theta_deg[0],theta_deg[1],phi_deg[0],phi_deg[1],this_weight))
		ax.grid()
		cbar=pylab.colorbar(imgplot)
		if fluxflag:
			cbar.set_label(r"n p$^{-1}$ cm$^{-2}$")
		else:
			cbar.set_label(r"n p$^{-1}$")
		#ax.plot(zap_x1,zap_y,color=[0.5,0.5,0.5],linewidth=3,linestyle='--')
		#ax.plot(zap_x2,zap_y,color=[0.5,0.5,0.5],linewidth=3,linestyle='--')
		ax.set_xlim([x_bins[0],x_bins[-1]])
		ax.set_ylim([y_bins[0],y_bins[-1]])
		#
		# 10
#		#
#		if   theta_bin ==0 and E_bin == 0:
#			cbar.set_clim(0, 1.5e-5) #5e-6)
#		elif theta_bin ==0 and E_bin == 1:
#			cbar.set_clim(0, 7.2e-6)
		#
		# 90
		#
#		if   theta_bin ==0 and E_bin == 0:
#			cbar.set_clim(2.1e-4, 4.1e-4) #5e-6)
#		elif theta_bin ==0 and E_bin == 1:
#			cbar.set_clim(0, 2.6e-4)
#		elif theta_bin ==0 and E_bin == 2:
#			cbar.set_clim(0, 2e-6)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		f.savefig('dist_e%d_theta%d'%(E_bin,theta_bin))
		pylab.show()


### spectrum plots 
f3 = plt.figure()
gs = gridspec.GridSpec(3, 4)
ax3 = f3.add_subplot(gs[0,:-1])
ax4 = f3.add_subplot(gs[1,:-1])
#ax5 = f3.add_subplot(gs[2,:-1])
ax6 = f3.add_subplot(gs[2,:-1])

cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=len(theta_bins)-1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

hist=[]
err=[]
area=(y_bins[-1]-y_bins[0])*(x_bins[-1]-x_bins[0])
for i in range(0,len(theta_bins)-1):
	sa = 1 #numpy.pi * ( theta_bins[i+1]*theta_bins[i+1] - theta_bins[i]*theta_bins[i] )
	if fluxflag:
		h = histograms_flux[i].values
	else:
		h = histograms_curr[i].values

	spec = numpy.divide(h,numpy.diff(histograms_flux[i].bins))
	colorVal = scalarMap.to_rgba(i)
	ax3.semilogx(histograms_curr[i].bins[:-1],numpy.divide(spec,sa),color=colorVal,label=r'$\theta$ = %4.2f - %4.2f'%(theta_bins[i]*180.0/numpy.pi,theta_bins[i+1]*180.0/numpy.pi),drawstyle='steps-post',linewidth=1)
	ax4.semilogx(histograms_curr[i].bins[:-1],numpy.divide(h,   sa),color=colorVal,label=r'$\theta$ = %4.2f - %4.2f'%(theta_bins[i]*180.0/numpy.pi,theta_bins[i+1]*180.0/numpy.pi),drawstyle='steps-post',linewidth=1)

handles, labels = ax3.get_legend_handles_labels()
ax3.legend(handles,labels,loc=1,prop={'size':12}, ncol=2, bbox_to_anchor=(1.4, 1.1))
#handles, labels = ax4.get_legend_handles_labels()
#ax4.legend(handles,labels,loc=1,prop={'size':12}, ncol=2, bbox_to_anchor=(1.4, 1.1))


f=open('%d.theta_spec'%this_sc,'w')
cPickle.dump([finebins,hist,err,y_bins,x_bins,theta_bins],f)
f.close()


if fluxflag:
	ax3.set_ylabel(r'$\Phi$(E) dE (n p$^{-1}$ cm$^{-2}$ MeV$^{-1}$ Str$^{-1}$)')
	ax4.set_ylabel(r'$\Phi$(E) (n p$^{-1}$ cm$^{-2} Str$^{-1}$$)')
	ax6.set_ylabel(r'$\Phi$(E) dE (n p$^{-1}$ cm$^{-2}$ MeV$^{-1}$ Str$^{-1}$)')
else:
	ax3.set_ylabel(r'$\Phi$(E) dE (n p$^{-1}$ MeV$^{-1}$ )')
	ax4.set_ylabel(r'$\Phi$(E) (n p$^{-1}$               )')
	ax6.set_ylabel(r'$\Phi$(E) dE (n p$^{-1}$ MeV$^{-1}$ )')
#ax5.set_ylabel(r'Cumulative Probability')
#ax5.set_ylim([0,1])
ax3.grid()
ax4.grid()
#ax5.grid()
ax6.grid()

ax3.set_xlim([1e-10,1e-4])
ax4.set_xlim([1e-10,1e-4])
ax6.set_xlim([1e-10,1e-4])

pylab.show()



import scipy.io
data_to_save={}
for i in range(0,len(theta_bins)-1):
	string = 'dist_%3.2E_%3.2E'%(theta_bins[i],theta_bins[i+1])
	data_to_save[string]=dist[0][i][phi_bin][:][:]
data_to_save['y_bins']=y_bins
data_to_save['x_bins']=x_bins
scipy.io.savemat('FOCUS_divergence.mat',data_to_save)
