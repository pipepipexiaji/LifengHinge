## usage: makePsf origami-90base.pdb  villin-noWater

puts "Usage: makePsf {inPdb outFile}"

proc makePsf {inPdb outFile} {
    ## initialize psfgen
    package require psfgen             ;# load psfgen
    resetpsf                           ;# Destroys any previous attempts
    psfcontext reset                   ;# Tosses out any previously declared topology files

    psfalias                           ;# tell psfgen what resnames to expect
    topology /Users/admin/Desktop/basic-info/top_all36_na.rtf ;# tells psfgen how atoms in a residue should be bound
    topology /Users/admin/Desktop/basic-info/toppar_water_ions_jejoong.str
    
    ## load inPdb
    set ID [mol new $inPdb]

    ## add atom selections to psfgen
    
    buildDNASegmentFromSel [atomselect $ID "nucleic"] SCF
    
    #buildDNASegmentFromSel [atomselect $ID "nucleic and not (segname SC0 SC1 SC2 SC3)"] STP
    buildDNASegmentFromSel [atomselect $ID "segname MGHH"] ION 
    #buildDNASegmentFromSel [atomselect $ID "resname TIP3"] WT
    #buildDNASegmentFromSel [atomselect $ID "name POT"] POT
    #buildDNASegmentFromSel [atomselect $ID "name CLA"] CLA

    ## psfgen commands
    guesscoord ;# guess the coordinates of missing atoms
    regenerate angles dihedrals ;# fixes problems with patching

    ## write psf and pdb that were added in psfgen's 
    writepdb $outFile.pdb
    writepsf $outFile.psf
    
    ## cleanup in VMD
    mol delete $ID
}


##########################################################################
#### procs used by makePsf (kindly ignore, unless you are interested) ####
##########################################################################
proc buildDNASegmentFromSel { selection {segSuffix PRO} } {
    ## add $sel to psfgen
    set tmpDir /var/tmp ;# somewhere to drop temporary files
    set ID [$selection molid] ;# current VMD mol
    

    ## psfgen works by loading a pdb for each segment (a unique 4-character identifier for a group of atoms that are usually bonded together e.g. a single component from a PDB)
    ##   these can usually be identified from the "chain" field of a pdb

    ## add each chain in $sel to psfgen as its own segment
    set segnames [lsort -unique [$selection get segname]]
    foreach segname $segnames {
	set seg ${segname}
	set sel [atomselect $ID "[$selection text] and segname $segname"]
	

	## write out temporary pdb for psfgen to read
	set tmpPdb $tmpDir/tmp.pdb
	$sel writepdb $tmpPdb
	
	## add segment to psfgen using psfgen commands
	
	if {$segSuffix == "SCF" } {
		segment $seg {
	    		first 5TER
	    		last 3TER
	    		pdb $tmpPdb
		}
		
		## Now we patch the RNA molecule made by psfgen to make DNA
		##     By default, psfgen makes DNA.
	
		foreach {patchSeg patchResId patchResName} [join [lsort -unique [$sel get {segname resid  resname}]]] {
			if {$patchResId == 1} {
            			patch DEO5 $patchSeg:$patchResId 
        		} else {
				patch DEOX $patchSeg:$patchResId
			}
            	
		}

		# But we also need to make the DNA periodic
        
		#set first [lindex [$sel get resid] 0]
		#set last [lindex [$sel get resid] end]
		#patch LKNA $seg:$last $seg:$first
		
		coordpdb $tmpPdb
	
	}

	if {$segSuffix == "STP" } {

		segment $seg {
	    		first 5TER
	    		last 3TER
	    		pdb $tmpPdb
		}
		
		## Now we patch the RNA molecule made by psfgen to make DNA
		##     By default, psfgen makes DNA.
	
		foreach {patchSeg patchResId patchResName} [join [lsort -unique [$sel get {segname resid  resname}]]] {
			if {$patchResId == 1} {
            			patch DEO5 $patchSeg:$patchResId 
        		} else {
				patch DEOX $patchSeg:$patchResId
			}
            	
		}	
		#patch LKNA $seg:28 $seg:29
		coordpdb $tmpPdb
	
	}

	if {$segSuffix != "SCF" && $segSuffix != "STP"} {

		segment $seg {
	    		first none
	    		last none
	    		pdb $tmpPdb
		}
		coordpdb $tmpPdb
	}

	if {$segSuffix == "DNA"} {
		## Now we patch the RNA molecule made by psfgen to make DNA
		##     By default, psfgen makes DNA.
	
		foreach {patchSeg patchResId patchResName} [join [lsort -unique [$sel get {segname resid  resname}]]] {

	    		

	    		patch DEOX $patchSeg:$patchResId
            	

		}

		# But we also need to make the DNA periodic
        
		set first [lindex [$sel get resid] 0]
		set last [lindex [$sel get resid] end]
		patch LKNA $seg:$last $seg:$first
	}

		## That was fun and easy

    }
    
    if {0} {
	
	patch LKNA P0:28 P0:29
	patch LKNA P1:28 P1:29
	patch LKNA P2:28 P2:29
	patch LKNA P3:28 P3:29
	patch LKNA P4:28 P4:29
	patch LKNA P5:28 P5:29
	patch LKNA P6:28 P6:29
	patch LKNA P7:28 P7:29

    }

}    

proc psfalias {} { ##
    ## Define common aliases for psfgen (so that it knows what residue and atom names mean)

    # Here's for nucleics
    foreach suff {"" 3 5} {
        pdbalias residue DA$suff ADE
        pdbalias residue DT$suff THY
        pdbalias residue DC$suff CYT
        pdbalias residue DG$suff GUA
        pdbalias residue A$suff ADE
        pdbalias residue T$suff THY
        pdbalias residue C$suff CYT
        pdbalias residue G$suff GUA
    }


    foreach bp { GUA CYT ADE THY URA } {
	pdbalias atom $bp "O5\*" O5'
	pdbalias atom $bp "C5\*" C5'
	pdbalias atom $bp "O4\*" O4'
	pdbalias atom $bp "C4\*" C4'
	pdbalias atom $bp "C3\*" C3'
	pdbalias atom $bp "O3\*" O3'
	pdbalias atom $bp "C2\*" C2'
	pdbalias atom $bp "O2\*" O2'
	pdbalias atom $bp "C1\*" C1'
	pdbalias atom $bp  OP1   O1P
	pdbalias atom $bp  OP2   O2P
    }
    pdbalias atom THY C7 C5M

    pdbalias atom ILE CD1 CD
    pdbalias atom SER HG HG1
    pdbalias residue HIS HSE

    # Heme aliases
    pdbalias residue HEM HEME
    pdbalias atom HEME "N A" NA
    pdbalias atom HEME "N B" NB
    pdbalias atom HEME "N C" NC
    pdbalias atom HEME "N D" ND

    # Water aliases
    pdbalias residue HOH TIP3
    pdbalias atom TIP3 O OH2
    pdbalias atom TIP3 OW OH2

    # Ion aliases
    pdbalias residue K POT
    pdbalias atom POT K POT
    pdbalias residue ICL CLA
    pdbalias atom ICL CL CLA
    pdbalias residue NA SOD
    pdbalias atom SOD NA SOD
    pdbalias residue CA CAL
    pdbalias atom CAL CA CAL 
    
    pdbalias atom ATP C1* C1'
    pdbalias atom ATP C2* C2'
    pdbalias atom ATP C3* C3'
    pdbalias atom ATP C4* C4'
    pdbalias atom ATP C5* C5'
    pdbalias atom ATP O2* O2'
    pdbalias atom ATP O3* O3'
    pdbalias atom ATP O4* O4'
    pdbalias atom ATP O5* O5'
}
