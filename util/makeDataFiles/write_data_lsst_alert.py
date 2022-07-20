# RKnop 2022-06-30 : Rewritten.  (I found the logic of the original a
#          bit hard to follow, so I wrote it afresh so that it can
#          straightforwardly keep track of instance variables, make
#          fresh alert data, and keep arrays of previous alert data.)
#
#          What's different about the output: populates the
#          prvDiaForcedSources array.  Now includes *all* previous
#          *detections* in prvDiaSources.  (Before, it put what should
#          have been in prvDiaForcedSources instead into prvDiaSources.)
#          Also, read the schema and build the base objects to go into
#          the alert directly from the .avsc file, not depending on any
#          example alerts in .json files.  Finally, alerts are now
#          written out without embedded schema; since we don't need to
#          stream to zads (as they will be *sent* schemaless), there's
#          no need to write the schema with each alert.

# Previous notes:
# Created Oct 22, 2021
# write data in lsst-alert format for broker test.
# [R.Kessler, G. Narayan, R.Hlozek, ]
#
# Jan 14 2022 G.Narayan - fix diaObject bug found by Rob K.
# Jan 15 2022 R.Kessler - minor cleanup; start working on reducing output
# Jan 31 2022 RK^2 fix bug setting alert_first_detect
# Feb 22 2022 RK - integreate ZPHOT_Q and check for 2nd hostgal
# Mar 30 2022 RK - fix setting alertId (not alertID) for all epochs
#                   (not just for FIRST_OBS)
# Mar 20 2022 RK - create & write simVersion string
# 
# Apr 09 2022 RK - few updates to run x10 faster:
#   + for sunset-MJD, replace astroplan call with grid search on
#       pre-computed list of sunset-MJDs at CTIO (x5 faster)
#   + use os.mkdir to create nite-dir instead os is.system (x2 faster)
#
# May 12 2022 RK - hostgal_z[_err] -> hostgal_zspec[_err] and add
#                  hostgal_zphot[_err]
#
# Jun 15 2022 RK - add HOSTKEY_RA[DEC] to list


import os, sys, copy, re, gzip, math, pathlib
import datetime, logging
import numpy

import fastavro
import fastavro.utils

import makeDataFiles_params as gpar
import makeDataFiles_util as util

_logger = logging.getLogger(__name__)
if not _logger.hasHandlers():
    _logout = logging.StreamHandler( sys.stderr )
    _logger.addHandler( _logout )
    _formatter = logging.Formatter( f'[%(asctime)s - %(levelname)s] - %(message)s',
                                    datefmt='%Y-%m-%d %H:%M:%S' )
    _logout.setFormatter( _formatter )
    _logger.propagate = False
    _logger.setLevel( logging.INFO )
    # _logger.setLevel( logging.DEBUG )
    

class LsstAlertWriter:

    PHOTFLAG_DETECT = 4096            # Read this from global data header?
    GZIP_ALERTS = True
    
    def __init__( self, args, config_data,
                  max_alerts_per_obj=2000, print_update_every=100, zp=31.4,
                  alert_day_name="NITE", lsst_site_name="CTIO", logger=_logger ):
        """Create an object that can write LSST alerts

        args — passed command-line args
        config_data — Something outside is expecting a couple of fields (n_alert_write, n_event_write) 
                       to be updated
        max_alerts_per_obj — Maximum alerts there will ever be for an object.  This is used to
                              generate unique SourceIDs, so make it big enough!  (Default 2000)
        print_every — Print status update to the console every time this many alerts are generated
                      (Default 100)
        zp — The assumed LSST zeropoint for reporting calibrated flux (Default 31.4)
        alert_day_name — (Default "NITE")
        lsst_site_name — (Default "CTIO") 
        """

        self.logger = logger
        self.t_start = None
        self.tot_n_alerts_written = 0
        self.next_print = print_update_every
        self.tot_n_events_processed = 0
        self.print_update_every = print_update_every
        self.max_alerts_per_obj = max_alerts_per_obj
        self.zp = zp
        self.scale_fluxcal = math.pow( 10. , 0.4 * (self.zp - gpar.SNANA_ZP) )
        self.keyname_substring_fluxcal = 'FLUXCAL'
        self.outdir = args.outdir_lsst_alert
        self.truthfile = args.outfile_alert_truth
        if self.truthfile:
            self.truth_dict = {
                'outfile':      pathlib.Path( self.truthfile ),
                'diaSourceId':  [],
                'snid':         [],
                'mjd':          [],
                'detect':       [],
                'true_gentype': [],
                'true_genmag':  []
            }
        else:
            self.truth_dict = None

        self.alert_day_name = alert_day_name
        self.lsst_site_name = lsst_site_name
        if args.mjd_sunset_file:
            self.mjd_sunset_dict = { 'mjd_file': args.mjd_sunset_file }
        else:
            self.mjd_sunset_dict = {}
            
        # Notice we're not making a copy of these; so, if the outside
        #   data changes, the data inside the structure here will change too
        self.args = args           
        self.config_data = config_data

        # Build some maps for converting things from SNANA variables to alert schema

        self.diaSource_map = {
            gpar.DATAKEY_SNID     : 'diaObjectId',
            gpar.DATAKEY_RA       : 'ra',
            gpar.DATAKEY_DEC      : 'decl',
            gpar.DATAKEY_NOBS     : 'nobs',       # in phot_raw, not header
            'MJD'                 : 'midPointTai',
            'BAND'                : 'filterName',
            'FLUXCAL'             : 'psFlux',
            'FLUXCALERR'          : 'psFluxErr'
        }

        self.diaForcedSource_map = {
            gpar.DATAKEY_SNID            : 'diaObjectId',
            gpar.DATAKEY_NOBS            : 'nobs',       # in phot_raw, not header
            'MJD'                        : 'midPointTai',
            'BAND'                       : 'filterName',
            'FLUXCAL'                    : 'psFlux',
            'FLUXCALERR'                 : 'psFluxErr',
        }

        # NOTE: this next map isn't complete.  However, gpar doesn't
        #  currently have everything that needs to be in the map.  We
        #  have to rebuild the map dynamically, as gpar will be updated
        #  as necessary; that will be done before it's used in
        #  append_hostgal_diaObject_map
        self.diaObject_map_base = {
            gpar.DATAKEY_SNID            : 'diaObjectId',
            gpar.DATAKEY_RA              : 'ra',
            gpar.DATAKEY_DEC             : 'decl',
            gpar.DATAKEY_MWEBV           : 'mwebv',
            gpar.DATAKEY_MWEBV_ERR       : 'mwebv_err',
            gpar.DATAKEY_zHEL            : 'z_final' ,
            gpar.DATAKEY_zHEL_ERR        : 'z_final_err',
            gpar.HOSTKEY_ELLIP           : 'hostgal_ellipticity',
            gpar.HOSTKEY_SQRADIUS        : 'hostgal_sqradius',
            gpar.HOSTKEY_SPECZ           : 'hostgal_zspec',
            gpar.HOSTKEY_SPECZ_ERR       : 'hostgal_zspec_err',
            gpar.HOSTKEY_PHOTOZ          : 'hostgal_zphot',
            gpar.HOSTKEY_PHOTOZ_ERR      : 'hostgal_zphot_err',
            gpar.HOSTKEY_RA              : 'hostgal_ra',
            gpar.HOSTKEY_DEC             : 'hostgal_dec',
            gpar.HOSTKEY_SNSEP           : 'hostgal_snsep',
        }

        # Read the schema
        # The namespace is elasticc.v0_9 (at least as of this writing)
        # The base schema is alert, but we also have embedded schema
        # diaSource, diaForcedSource, diaNondetectionLimit, and diaObject
        self.alert_schema = fastavro.schema.parse_schema( fastavro.schema.load_schema( args.lsst_alert_schema ) )
        # There may be a nice way to get the namespace with fastavro,
        # but I don't know it, so I'm just going to parse it out of the name
        match = re.search( '^(.*)\.alert$', self.alert_schema['name']  )
        if match is None:
            raise RuntimeError( f"Failed to parse schema name {schema['name']} for (.*)\.alert" )
        namespace = match.group(1)
        self.object_schema = self.alert_schema['__named_schemas'][f'{namespace}.diaObject']
        self.source_schema = self.alert_schema['__named_schemas'][f'{namespace}.diaSource']
        self.forced_source_schema = self.alert_schema['__named_schemas'][f'{namespace}.diaForcedSource']

        # Build the simVersion string
        vals = util.extract_sim_readme_info( args.snana_folder,
                                             ['TIME_START',
                                              'SNANA_VERSION',
                                              'SIMLIB_FILE'] )
        cadence = os.path.basename( vals['SIMLIB_FILE'] ).rsplit('.',1)[0]
        t_start = str(vals['TIME_START']).split()[0]          # just y/m/d

        self.simVersion = ( f"date({t_start})_"
                            f"snana({vals['SNANA_VERSION']})_"
                            f"cadence({cadence})_"
                            f"schema({namespace})" )

        self.logger.info( f"Initialized LsstAlertWriter for simVersion {self.simVersion}" )


    def append_hostgal_diaObject_map( self ):
        """Append dynamic host galaxy columns.

        gpar will have been updated in real time outside of this class,
        so the diaObject map needs to be updated.
        """

        diaObject_map = self.diaObject_map_base.copy()
        
        for prefix in [ 'HOSTGAL_MAG', 'HOSTGAL_MAGERR' ]:
            for band in list( gpar.SURVEY_INFO['FILTERS']['LSST'] ):
                diaObject_map[ f"{prefix}_{band}" ] = f"{prefix.lower()}_{band}"

        for key in gpar.DATAKEY_LIST_ZPHOT_Q:
            if "HOSTGAL2" in key: continue
            diaObject_map[ key ] = key.lower()

        # For each HOSTGAL_XXX key, add another key with HOSTGAL2_XXX
        hostkeys = [ key for key in diaObject_map if gpar.HOSTKEY_BASE in key ]
        for key in hostkeys:
            newkey = util.key_hostgal_nbr( key, 2 )
            newval = util.key_hostgal_nbr( diaObject_map[key], 2 )
            diaObject_map[ newkey ] = newval

        return diaObject_map
        

    def get_event_dict_value( self, data_event_dict, varname, obsnum=None ):
        """Tries to pull the value of varname out of the right subdictionary of data_event_dict"""

        head_raw = data_event_dict['head_raw']
        head_calc = data_event_dict['head_calc']
        phot_raw = data_event_dict['phot_raw']
        value = None

        if varname in head_raw:
            value = int(head_raw[varname]) if varname == gpar.DATAKEY_SNID else head_raw[varname]
        elif varname in head_calc:
            value = head_calc[varname]
        elif varname in phot_raw:
            if ( isinstance( phot_raw[varname], numpy.ndarray ) or
                 isinstance( phot_raw[varname], list ) ):
                value = phot_raw[varname][obsnum]
            else:
                value = phot_raw[varname]
        else:
            raise ValueError( f"Didn't find {varname} in data_event_dict" )

        return value
            
        
    def make_object_alert( self, data_event_dict ):
        diaObject_map = self.append_hostgal_diaObject_map()
        diaObjectAlert = fastavro.utils.generate_one( self.object_schema )
        head_raw = data_event_dict['head_raw']
        for data_key in diaObject_map:
            diaObjectAlert[ diaObject_map[data_key] ] = \
                self.get_event_dict_value( data_event_dict, data_key )

        return diaObjectAlert

    def make_source_alert( self, data_event_dict, obsnum, diaSourceId ):
        diaSourceAlert = fastavro.utils.generate_one( self.source_schema )
        head_raw = data_event_dict['head_raw']
        for data_key in self.diaSource_map:
            diaSourceAlert[ self.diaSource_map[data_key] ] = \
                self.get_event_dict_value( data_event_dict, data_key, obsnum )

        diaSourceAlert['snr'] = diaSourceAlert['psFlux'] / diaSourceAlert['psFluxErr']
        # I don't think there's anything in SNANA that corresponds to ccdFVisit
        diaSourceAlert['ccdVisitId'] = -1
        diaSourceAlert['diaSourceId'] = diaSourceId
            
        return diaSourceAlert

    def make_forced_source_alert( self, data_event_dict, obsnum, diaForcedSourceId ):
        diaForcedSourceAlert = fastavro.utils.generate_one( self.forced_source_schema )
        head_raw = data_event_dict['head_raw']
        for data_key in self.diaForcedSource_map:
            diaForcedSourceAlert[ self.diaForcedSource_map[data_key] ] = \
                self.get_event_dict_value( data_event_dict, data_key, obsnum )

        # I don't know what in SNANA should correspond to totFlux
        diaForcedSourceAlert[ 'totFlux' ] = -1.
        diaForcedSourceAlert[ 'totFluxErr' ] = 1.
        diaForcedSourceAlert[ 'diaForcedSourceId' ] = diaForcedSourceId

        return diaForcedSourceAlert
        
    def outdir_for_nite( self, mjd ):
        nite = int( util.get_sunset_mjd( mjd, self.lsst_site_name,
                                         self.mjd_sunset_dict ) )
        outdir = pathlib.Path( self.outdir ) / f"{self.alert_day_name}{nite}"
        if not outdir.exists():
            outdir.mkdir( mode=0o770, parents=True, exist_ok=True )
        return outdir
        
    def write_alerts_for_event( self, data_event_dict ):
        """Write out all alerts for a single event (i.e. object).

        data_event_dict — dictionary with all the snana information about all
                          the photometry etc. of this one event.

        Writes out an alert for all events that are flagged as detected.
        Includes all previous detected events in prvDiaSources.
        Includes some previous events (detected or not) in
        prvDiaForcedSources *if* the time of the detection is at least
        12 hours later than the time of the first detection; the range
        of time included in prvDiaForcedSources is limited.  (Usually:
        starting 30 days before the first detection.  But, there's a
        different mode if you set --peakmjd-range.)
        """

        if self.t_start is None:
            self.t_start = datetime.datetime.now()
        
        head_raw  = data_event_dict['head_raw']
        head_calc = data_event_dict['head_calc']
        head_sim  = data_event_dict['head_sim']
        phot_raw  = data_event_dict['phot_raw']

        # Adjust for zeropoint
        for key in self.diaSource_map:
            if self.keyname_substring_fluxcal in key:
                phot_raw[key] = [ f * self.scale_fluxcal for f in phot_raw[key] ]

        snid = int( head_raw[gpar.DATAKEY_SNID] )
        nobs = phot_raw[ gpar.DATAKEY_NOBS]
        true_gentype = head_sim[ gpar.SIMKEY_TYPE_INDEX ]
        
        # Make a list of dates, dettecctions, and magnitudes
        mjds = data_event_dict['phot_raw']['MJD']
        photflags = data_event_dict['phot_raw']['PHOTFLAG']
        detects = [ ( photflag & self.PHOTFLAG_DETECT ) > 0
                    for photflag in photflags ]
        true_genmags = data_event_dict['phot_raw'][ gpar.VARNAME_TRUEMAG ]

        # Figure out the first and last day to include forced photometry
        firstdetect_mjd = head_calc[ gpar.DATAKEY_MJD_DETECT_FIRST ]
        lastdetect_mjd = head_calc[ gpar.DATAKEY_MJD_DETECT_LAST ]
        if self.args.nite_detect_range:
            force_start = firstdetect_mjd - 30
            force_end = lastdetect_mjd
        elif self.args.peakmjd_range:
            mjd_ref = head_calc[ GPAR.DATAKEY_PEAKMJD ]
            force_start = mjd_ref - 50
            force_end = mjd_ref + 100
        else:
            force_start = -1e32
            force_end = 1e32

        self.logger.debug( f"Object {snid}, force_start={force_start:.2f}, "
                           f"force_end={force_end:.2f}, "
                           f"firstdetect_mjd={firstdetect_mjd:2f}, "
                           f"mjds[0]={mjds[0]}" )
        
        # The diaObject information will be the same for everything
        diaObject = self.make_object_alert( data_event_dict )
        
        prvSources = []
        prvForcedSources = []
        # Go through all the events, building up prvForcedSources and
        # generating alerts as necessary
        # NOTE : Assuming here that mjd is coming in order!
        # I hope that's right.  (It seems to be.)
        # tmpmjd = [ i for i in mjds ]
        # tmpmjd.sort()
        # tmpmjd = numpy.array( tmpmjd )
        # if not numpy.all( mjds == tmpmjd ):
        #     self.logger.warning( f"MJDS not sorted for diaObjectId = {diaObject.diaObjectid}" )
        # else:
        #     self.logger.info( f"MJDs are sorted for diaObjectId = {diaObject.diaObjectid}" )
        for obsnum in range( nobs ):
            diaSourceId = snid * self.max_alerts_per_obj + obsnum
            
            diaForcedSource = None
            if ( ( mjds[ obsnum ] >= force_start )
                 and
                 ( mjds[ obsnum ] <= force_end ) ):
                diaForcedSource = self.make_forced_source_alert(
                    data_event_dict, obsnum, diaSourceId
                )

            # Only write to the truth table if it's a forced source
            #  or it's detected; otherwise, it will never show
            #  up anywhere in any alert.
            if self.truth_dict is not None:
                if detects[obsnum] or ( diaForcedSource is not None ):
                    self.truth_dict['diaSourceId'].append( diaSourceId )
                    self.truth_dict['snid'].append( snid )
                    self.truth_dict['mjd'].append( mjds[obsnum] )
                    self.truth_dict['detect'].append( detects[obsnum] )
                    self.truth_dict['true_gentype'].append( true_gentype )
                    self.truth_dict['true_genmag'].append( true_genmags[obsnum] )
                
            diaSource = None
            if detects[ obsnum ]:
                diaSource = self.make_source_alert(
                    data_event_dict, obsnum, diaSourceId
                )
                
                alert = fastavro.utils.generate_one( self.alert_schema )
                alert['diaObject'] = diaObject
                alert['diaSource'] = diaSource
                alert['prvDiaSources'] = prvSources
                alert['alertId'] = diaSource['diaSourceId']    # Hey, it's unique
                
                # Require a 12-hour delay between detection and
                # inclusion of forced photometry
                if ( mjds[ obsnum ] - 12 > firstdetect_mjd ):
                    alert['prvDiaForcedSources'] = prvForcedSources

                daystr = f"mjd{mjds[obsnum]:.4f}"
                objstr = f"obj{snid}"
                srcstr = f"src{diaSourceId}"

                if self.GZIP_ALERTS:
                    gz = ".gz"
                    openfunc = lambda x: gzip.open( x, mode='wb', compresslevel=9 )
                else:
                    gz = ""
                    openfunc = lambda x: open( x, mode='wb' )
                outfile = ( self.outdir_for_nite( mjds[obsnum] ) /
                            f"alert_{daystr}_{objstr}_{srcstr}.avro{gz}" )
                with openfunc( outfile ) as ofp:
                    fastavro.write.schemaless_writer( ofp, self.alert_schema, alert )

                self.tot_n_alerts_written += 1
                # Something outside needs the number of alerts
                self.config_data[ 'n_alert_write' ] = self.tot_n_alerts_written
                prvSources.append( diaSource )

            if diaForcedSource is not None:
                prvForcedSources.append( diaForcedSource )

        self.tot_n_events_processed += 1
        # Something outside needs the number of events
        self.config_data[ 'n_event_write' ] = self.tot_n_events_processed

        # Log every so often
        
        if self.tot_n_alerts_written >= self.next_print:
            self.print_writerate()
            while self.next_print <= self.tot_n_alerts_written:
                self.next_print += self.print_update_every


    def print_writerate( self ):
        deltat = datetime.datetime.now() - self.t_start
        secs = deltat.total_seconds()
        if secs > 7200:
            t = f"{secs/3600:.2f} hours"
        elif secs > 120:
            t = f"{secs/60.:.2f} minutes"
        else:
            t = f"{secs:.1f} seconds"
        rate = self.tot_n_alerts_written / secs
                
        self.logger.info( f"Wrote {self.tot_n_alerts_written} alerts "
                          f"for {self.tot_n_events_processed} events "
                          f"in {t} ({rate:.1f} alerts/second)" )

    def write_truth( self ):
        outfile = ( self.truth_dict['outfile'].parent /
                    f"{self.truth_dict['outfile'].name}.gz" )
        self.logger.info( f"Writing {len(self.truth_dict['diaSourceId'])} "
                          f"entries to truth table {outfile}" )
        with gzip.open( outfile, mode='wt' ) as ofp:
            ofp.write(f"SourceID, SNID, MJD, DETECT, "
                      f"TRUE_GENTYPE, TRUE_GENMAG\n")
            for obsnum in range(len(self.truth_dict['diaSourceId'])):
                ofp.write( f"{self.truth_dict['diaSourceId'][obsnum]}, "
                           f"{self.truth_dict['snid'][obsnum]}, "
                           f"{self.truth_dict['mjd'][obsnum]:.4f}, "
                           f"{1 if self.truth_dict['detect'][obsnum] else 0}, "
                           f"{self.truth_dict['true_gentype'][obsnum]}, "
                           f"{self.truth_dict['true_genmag'][obsnum]:.3f}\n" )
        
    def finalize( self ):
        self.print_writerate()
        self.write_truth()
        
        
