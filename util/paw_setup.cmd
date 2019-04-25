#!/bin/csh
#
# setup PAW macros to look at histograms
# produced by snana
# This script copies/updates your .pawlogon.kumac file,
# and it copies a few *.kumac files to ~/kumacs
# It only copies files that don't exist.
#
#
# Oct 20, 2006:  add 'opt zfl' and 'opt nbox' for make ps files
#
# Oct 12, 2007:  fix path bug : inputs/paw instead of just input
#
# May 14, 2008: add "FILECASE KEEP" to .pawlogon file
#
# Jun 2, 2008: copy .kumac files from $SNANA_DIR/kumacs
#              instead of inputs/paw
#
# ----------------------------------
  
# first create ~/kumacs directory

  cd
  if ( -d  kumacs ) then
     echo " ~/kumacs already exists. "
  else
    echo " Creating kumacs directory "
    mkdir kumacs 
  endif

  cd kumacs
  pwd

  set kumaclist = ( "hpl.kumac"  "overall.kumac"  "snana.kumac" "journal.kumac" )

  foreach kumac ( $kumaclist )
    if ( -e $kumac  ) then
      echo "    $kumac already exists. "
    else
      cp $SNANA_DIR/kumacs/$kumac .
      echo "    $kumac has been copied. "
    endif
  end

  echo ' '

# Now add stuff to .pawlogon so that hpl works

  cd
  set tmpfile = .pawlogon.kumac
  if ( -e $tmpfile ) then
     echo "  $tmpfile already exists. Removing it."
     rm $tmpfile
  endif
     echo "  Creating $tmpfile "
     touch $tmpfile

     echo "FILECASE KEEP" >> $tmpfile
     echo "set mtyp 20" >> $tmpfile

  set bla = `grep "opt zfl" $tmpfile`
  if ( $#bla == 0 ) then
     echo "opt zfl" >> $tmpfile
  endif
  set bla = `grep "opt nbox" $tmpfile`
  if ( $#bla == 0 ) then
     echo "opt nbox" >> $tmpfile
  endif


  set bla = `grep "auto" $tmpfile`
  if ( $#bla == 0 ) then
     echo " "  >> $tmpfile
     echo "* ----------------------------------- "  >> $tmpfile
     echo "macro/default '.,~/kumacs' -auto"  >> $tmpfile
  endif


  set bla = `grep "exec hpl" $tmpfile`
  if ( $#bla == 0 ) then
     echo "alias/create hpl 'exec hpl' c" >> $tmpfile
  endif

  set bla = `grep "exec hcopy" $tmpfile`
  if ( $#bla == 0 ) then
    echo "alias/create hcopy 'exec hpl#hcopy' c"  >> $tmpfile
  endif

  set bla = `grep "exec hfit" $tmpfile`
  if ( $#bla == 0 ) then
     echo "alias/create hfit 'exec hpl#hfit' c"   >> $tmpfile
  endif

  set bla = `grep "translation" $tmpfile`
  if ( $#bla == 0 ) then
    echo "alias/translation on"    >> $tmpfile
  endif

##########
  exit
##########


