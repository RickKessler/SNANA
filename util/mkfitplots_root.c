#include <iostream.h>
#include <TString.h>
#include <TMath.h>
int mkfitplots_root(char infile[200]="myfile.root", int SNid=-1, char options[100]="defaults", float ylimit1=20., float ylimit2=26.)
{
//
//  Eve Kovacs -- May 2013
//
//SNid=-1 - plots first SN in file
//SNid=1054 - plots SN with cid=1054 (or any other positive cid)
//SNid=-10  - plots first 10 (or any other negative number SNe in file)
//SNid=0 - plots last SN in file

//Plot options are comma separated and concatenated in any order (eg, "mutiscale,fit,print1" )
//"defaults" - row format, plot data and best-fit, auto-scale to same scale, plot flux(+ve & -ve), no diagnostic prints, no cuts, one pdf
//"snana" - defaults, but with pdf filename = infile.pdf (no annotations for option choices)
//"pdf" - output pdf file (default)
//"onepdf" - all plots written to same pdf file (default)
//"multipdf" - separate pdfs for each SN
//"gif" - output gif instead of pdf
//"onegif" - all plots written to same gif file (default for gif option)
//"multigif" - separate gifs for each SN
//"png" - output png instead of pdf
//"onepng" - all plots written to same png file (default for png option)
//"multipng" - separate pngs for each SN
//"multiscale" - separate scales for each filter plot
//"samescale" - same scale for each filter plot for a given SN (default)
//"onescale" - one scale for all filter plots
//"flux" - plot flux on y-axis (default)
//"posflux" - plot +ve fluxes only (ignore -ve fluxes)
//"redline" - plot a red line at y=0
//"mag" - plot magnitude on y-axis
//"auto" - auto-scale plots (default)
//"minmax" or "maxmin" - force scale on all plots to be between supplied values of ylimit1 and ylimit2 (flux or (-ve) magnitude,
//"min" - force minimum scale on all plots to be supplied value of  ylimit1 (flux or (-ve) magnitude, as requested)
//"max" - force maximum scale on all plots to be supplied value of  ylimit1 (flux or (-ve) magnitude, as requested)
//"row" - plots in horizontal format (default)
//"column" - plots in vertical format
//"data" - plot only data
//"fit" - plot only best fit
//"print1" - turn on some diagnostic printouts
//"print2" - turn on more diagnostic printouts
//"print3" - turn on all diagnostic printouts
//"SNR_a_b_c_d" - applies SNR cuts a,b,c,d,.. (any number up to 10) to flux-data points; skips SNe that fail cuts  
//"merge" - plots for each filter are not separated on page (has cosmetic defects) 

//Notes:
//Plots of flux have default minimum of 0, even if auto-scaling is on
//Minmax values can be in any order 


gROOT->Reset();
// Good setup for 1x4 plot
Int_t textfont=62;  //10*fontnumber + precision
Float_t textsize=.06; 
gStyle->SetOptStat(1111100); gStyle->SetHistLineWidth(3); gStyle->SetFuncWidth(3); 
gStyle->SetMarkerStyle(22); gStyle->SetMarkerSize(1.0); gStyle->SetMarkerColor(1); 
gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetOptFit(111); 
gStyle->SetTitleTextColor(1);  gStyle->SetTitleColor(1,"X"); gStyle->SetTitleColor(1,"Y"); 
gStyle->SetTitleColor(1,"X"); gStyle->SetTitleColor(1,"Y"); gStyle->SetTitleColor(1,"Z");
gStyle->SetTitleSize(textsize,"X"); gStyle->SetTitleSize(textsize,"Y"); gStyle->SetTitleSize(textsize,"Z");
gStyle->SetTitleOffset(1.0,"X"); gStyle->SetTitleOffset(1.0,"Y"); 
gStyle->SetPadBottomMargin(0.2); 
gStyle->SetPadLeftMargin(0.17);
gStyle->SetPadTopMargin(2.0); 
gStyle->SetPalette(1,0);
gStyle->SetTextFont(textfont);
gROOT->ForceStyle();
gStyle->SetCanvasColor(0);
gStyle->SetPadColor(0);
gStyle->SetFrameBorderMode(0);
gStyle->SetFillColor(0);
TGaxis::SetMaxDigits(3);

 Int_t i,j,k,l,n;
 char *pchr, *pchr1, *pchr2, *pchr3;
 Int_t const nchar=250;
 char filetype[nchar];
 Float_t ylim_lo=0., ylim_hi=0., maglim_lo_snana=26., fluxlimit;

 cout<<""<<endl;
 cout<<"Selected plotting options are "<<options<<"\n"<<endl;

 char printoptions[nchar],layout[nchar],outputoption[nchar],scaleoption[nchar],yscale[nchar],ylimits[nchar];
 sprintf(layout,"row");
 sprintf(outputoption,"one");
 sprintf(scaleoption,"same");
 sprintf(yscale,"flux");
 sprintf(printoptions,"");
 sprintf(filetype,"pdf");
 sprintf(ylimits,"auto");

//Force turn on of other options if selected
 if(pchr=strstr(options,"print1")!= NULL) sprintf(printoptions,"print1"); 
 if(pchr=strstr(options,"print2")!= NULL) sprintf(printoptions,"print1,print2");
 if(pchr=strstr(options,"print3")!= NULL) sprintf(printoptions,"print1,print2,print3");
 if(pchr=strstr(options,"debug")!= NULL) sprintf(printoptions,"debug,print1,print2,print3");
 if(pchr=strstr(options,"row")!= NULL) sprintf(layout,"row");
 if(pchr=strstr(options,"column")!= NULL) sprintf(layout,"column");
 if(pchr=strstr(options,"pdf")!= NULL) sprintf(filetype,"pdf");
 if(pchr=strstr(options,"gif")!= NULL) sprintf(filetype,"gif");
 if(pchr=strstr(options,"png")!= NULL) sprintf(filetype,"png");
 //if(pchr=strstr(options,"svg")!= NULL) sprintf(filetype,"svg");
 if(pchr=strstr(options,"onepdf")!= NULL) sprintf(outputoption,"one");
 if(pchr=strstr(options,"multipdf")!= NULL) sprintf(outputoption,"multi");
 if(pchr=strstr(options,"onegif")!= NULL) sprintf(outputoption,"one");
 if(pchr=strstr(options,"multigif")!= NULL) sprintf(outputoption,"multi");
 if(pchr=strstr(options,"onepng")!= NULL) sprintf(outputoption,"one");
 if(pchr=strstr(options,"multipng")!= NULL) sprintf(outputoption,"multi");
 //if(pchr=strstr(options,"onesvg")!= NULL) sprintf(outputoption,"one");
 //if(pchr=strstr(options,"multisvg")!= NULL) sprintf(outputoption,"multi");
 if(pchr=strstr(options,"multiscale")!= NULL) sprintf(scaleoption,"multi");
 if(pchr=strstr(options,"samescale")!= NULL) sprintf(scaleoption,"same");
 if(pchr=strstr(options,"onescale")!= NULL) sprintf(scaleoption,"one");
 if(pchr=strstr(options,"flux")!= NULL) sprintf(yscale,"flux");
 if(pchr=strstr(options,"mag")!= NULL) sprintf(yscale,"mag");
 if(pchr=strstr(options,"auto")!= NULL) sprintf(ylimits,"auto");
 if(pchr=strstr(options,"posflux")!= NULL){
   fluxlimit=0.;            //suppress negative flux
 }
 else {
   fluxlimit=-99999.;
 }
 if(strstr(options,"minmax")!= NULL || strstr(options,"maxmin")!= NULL) {  //force limits
   sprintf(ylimits,"minmax");
   ylim_lo=TMath::Min(ylimit1,ylimit2);
   ylim_hi=TMath::Max(ylimit1,ylimit2);
 }
 elseif(pchr1=strstr(options,"max")!= NULL) {
   sprintf(ylimits,"maxonly");
   ylim_hi=ylimit1; //for flux
   ylim_lo=ylimit1; //for mag
 }
 elseif(pchr1=strstr(options,"min")!= NULL) {
   sprintf(ylimits,"minonly");
   ylim_lo=ylimit1; //for flux
   ylim_hi=ylimit1; //for mag 
 }
 if(pchr=strstr(options,"snana")!= NULL) {   //override other selections if snana selected
   /*
    if(pchr=strstr(yscale,"mag")!= NULL) {
      if(strstr(ylimits,"auto")!= NULL) {
        sprintf(ylimits,"minonly");
        ylim_hi=maglim_lo_snana;
      }
    }
   */
    sprintf(outputoption,"onesnana");
 }
 if(pchr=strstr(options,"defaults")!= NULL) {
   sprintf(layout,"row");
   sprintf(outputoption,"one");
   sprintf(scaleoption,"same");
   sprintf(yscale,"flux");
   sprintf(printoptions,""); 
   sprintf(ylimits,"auto");
 }

 Int_t NSNe=1;  char SNcid[nchar]='';
 if (SNid >= 0) {
   sprintf(SNcid,"%d",SNid);
 }
 else {
   sprintf(SNcid,"-%d",abs(SNid));
   NSNe = TMath::Abs(SNid);
 }
 if (SNid <= 0) cout<<"Plotting lightcurves for "<<NSNe<<" SNe in each tree\n"<<endl;
 if (SNid > 0)  cout<<"Plotting lightcurves for SNcid "<<SNcid<<" in each tree\n"<<endl; 
 
 Int_t const Nsnrmax=10; Int_t SNRmax[Nsnrmax]={0}, isnr, lsnr, endsnr=0, Nsnr=0; 
 char snroption[nchar]=""; char tmp[nchar]="";
 if(pchr=strstr(options,"SNR_")!= NULL) {
   i=1;
   lsnr=(strlen(strstr(options,"SNR_"))-4);
   isnr=strlen(options)-lsnr-1;  //string index starts at 0
   while (i<=lsnr && endsnr==0 && Nsnr<Nsnrmax) {
     if(options[isnr+i] == ',') endsnr=1;
     if(pchr=strstr(printoptions,"debug")!= NULL) cout<<i<<","<<options[isnr+i]<<endl;
     if (options[isnr+i]!='_' && options[isnr+i] != ',') {
       sprintf(tmp,"%s%c",tmp,options[isnr+i]);
     }
     else {
       SNRmax[Nsnr]=atoi(tmp);
       if(pchr=strstr(printoptions,"debug")!= NULL) cout<<"Nsnr,SNRmax[Nsnr] "<<Nsnr<<","<<SNRmax[Nsnr]<<endl;
       Nsnr=Nsnr+1;
       sprintf(tmp,"");
     }
     i=i+1;
   }
 }
 if(pchr=strstr(printoptions,"debug")!= NULL) cout<<"Nsnr,tmp,endsnr "<<Nsnr<<","<<tmp<<","<<endsnr<<endl;
 if (strlen(tmp)>0) {
   SNRmax[Nsnr]=atoi(tmp); //save last value if  SNR option is at end of options string 
   Nsnr=Nsnr+1;
 }
 if (Nsnr>0) {
   sprintf(snroption,"SNR");
   cout<<"Using "<<Nsnr<<" SNR cuts :"<<endl;
      for (j=0;j<Nsnr;j++){
	cout<<" #"<<j+1<<" = "<<SNRmax[j]<<endl;
        sprintf(snroption,"%s_%d",snroption,SNRmax[j]);
      }
      cout<<" "<<endl;
 }

 char outputfile[nchar],openoutputfile[nchar],closeoutputfile[nchar];
 TFile *f = new TFile(infile,"READ");
 if(pchr=strstr(printoptions,"debug")!= NULL) f->ls();

 char survey[nchar]=""; char survey_filters[nchar]=""; char filterstring[nchar]="";
 Int_t const Nlcpakmax=10;
 Int_t nlcpak, temp; char tempstring[nchar];
 char snlcpak_treenames[nchar]=""; char snlcpak_tree[Nlcpakmax][nchar]={""}; 

 //TTree *global = GLOBAL;
 TTree *global = (TTree*)f->Get("GLOBAL");
 if(pchr=strstr(printoptions,"print3")!= NULL) global->Print();
 global->SetBranchAddress("SURVEY", survey);
 global->SetBranchAddress("SURVEY_FILTERS", survey_filters);
 global->SetBranchAddress("FILTERSTRING", filterstring);
 global->SetBranchAddress("NLCPAK",&nlcpak);
 global->SetBranchAddress("SNLCPAK_TREENAMES",snlcpak_treenames);

 Int_t nglobal = Int_t(global->GetEntriesFast());
 if(pchr=strstr(printoptions,"print3")!= NULL || nglobal > 1) cout<<"Number of entries in GLOBAL tree = "<<nglobal<<endl;
 global->GetEntry();
 cout<<" Survey filters for "<<survey<<" are: "<<survey_filters<<endl;
 cout<<" Filterstring is: "<<filterstring<<endl;
 cout<<" Number of SNLCPAK trees: "<<nlcpak<<endl;
 cout<<" SNLCPAK treenames: "<<snlcpak_treenames<<endl;  

//Setup Filter Ids and Pointer Arrays
 Int_t const Nfilt_max=80;
 Int_t Nfilters,Nfiltermax;
 char FilterId[Nfilt_max][nchar]={""};  //Filter Identifier array for use in titles etc
 Int_t IfilterId[Nfilt_max]={0}, Ifilt_min=Nfilt_max, Ifilt_max=0;

 Nfilters=strlen(survey_filters);
 Nfiltermax=strlen(filterstring);
 if (Nfiltermax > Nfilt_max) cout<<"Oops: Maximum number of filters exceeded!"<<endl;
 for (i=0;i<Nfilters;i++) {
   sprintf(FilterId[i],"%c",survey_filters[i]);   //write filters to 1-char string array
   if(pchr=strstr(printoptions,"debug")!= NULL) cout<<survey_filters[i]<<","<<FilterId[i]<<endl;
   if (pchr=strstr(filterstring,FilterId[i])!=NULL) {
     IfilterId[i]=Nfiltermax-strlen(strstr(filterstring,FilterId[i])); //Rick starts at 1
     if(pchr=strstr(printoptions,"print3")!= NULL) cout<<"Check of valid FilterID's: "<<FilterId[i]<<","<<IfilterId[i]<<endl;
     Ifilt_max=TMath::Max(Ifilt_max,IfilterId[i]);
     Ifilt_min=TMath::Min(Ifilt_min,IfilterId[i]);
   }
   else {
     cout<<"Oops: FilterId[i] = "<<FilterId[i]<<" not found in "<<filterstring<<endl;
   }
 }
 if(pchr=strstr(printoptions,"print1")!= NULL) cout<<"Minimum filter id = "<< Ifilt_min<<endl;
 if(pchr=strstr(printoptions,"print1")!= NULL) cout<<"Maximum filter id = "<< Ifilt_max<<"\n"<<endl;
 
//Remove leading blanks and other blanks to parse snlcpak_treenames
 if (pchr=strstr(snlcpak_treenames," ")!= NULL) {   //check for blanks
  Int_t newstring=0;
  n=0;
  for (i=0;i<strlen(snlcpak_treenames);i++) { //skip blanks
   if (snlcpak_treenames[i] != ' ') {               //if (snlcpak_treenames[i] != " ") doesn't work
     sprintf(snlcpak_tree[n],"%s%c",snlcpak_tree[n],snlcpak_treenames[i]);
     if(pchr=strstr(printoptions,"debug")!= NULL) cout<<"Check of parsing tree names"<<snlcpak_tree[n]<<","<<strlen(snlcpak_tree[n])<<endl;
     newstring = 1;
   }
   else {
     if (newstring ==1) { // first blank after previous string
       newstring=0;
       n++;
     }
   }
  }
  if (newstring ==1) n++; //no ending blanks found -- adjust count of parsed treenames
  if (n != nlcpak) cout<<"Oops: mismatch between snlcpak_treenames and nlcpak: # of treenames found  = "<<n<<endl; 
 }

//loop through all trees and selected supernovae

 Int_t const Nplotmax=2000;    //Check array sizes
 if(NSNe*nlcpak>Nplotmax) {
     NSNe=Nplotmax/nlcpak;
     cout<<"Warning: Maximum number ("<<Nplotmax<<") of plots exceeded!"<<endl; 
     cout<<"Only plotting lightcurves for first "<<NSNe<<" SNe\n"<<endl;
 }

 char Ccid[nchar];
 const Int_t Ntot_max=1000;  //Maxm number of observations
 Int_t Ntot;
 Float_t Fluxcal_data[Ntot_max],Fluxcal_data_error[Ntot_max],Tobs_data[Ntot_max],Fitchi2[Ntot_max];
 Float_t  Fluxcal_fit[Ntot_max],Fluxcal_fit_error[Ntot_max],Tobs_fit[Ntot_max];
 Int_t Ifiltobs[Ntot_max]; char Cfiltobs[Ntot_max][nchar];
 Int_t Iflagdata[Ntot_max];
 Int_t Reject[Ntot_max];
 Float_t Epoch[Ntot_max][Nfilt_max], Flux[Ntot_max][Nfilt_max], Flux_error[Ntot_max][Nfilt_max];
 Float_t Epochmin=999., Epochmax=-999., Fluxmax[Nplotmax][Nfilt_max]={0}, Flux_max[Nplotmax]={0.};
 Float_t DataEpochmin=999., DataEpochmax=-999., DataFluxmax[Nplotmax][Nfilt_max]={0}, DataFlux_max[Nplotmax]={0.};
 Float_t FitEpochmin=999., FitEpochmax=-999., FitFluxmax[Nplotmax][Nfilt_max]={0}, FitFlux_max[Nplotmax]={0.};
 Float_t FluxMAX=0., DataFluxMAX=0., FitFluxMAX=0.;
 Float_t Fluxmin[Nplotmax][Nfilt_max], Flux_min[Nplotmax];
 Float_t DataFluxmin[Nplotmax][Nfilt_max], DataFlux_min[Nplotmax];
 Float_t FitFluxmin[Nplotmax][Nfilt_max], FitFlux_min[Nplotmax];
 Float_t FluxMIN=99999., DataFluxMIN=99999., FitFluxMIN=99999.;
 Float_t zeropt=27.5;
 Double_t Band_peakmjd[Nfilt_max]={0};

 for (i=0;i<Nplotmax;i++){      //Setup minimum-value arrays
     for (j=0;j<Nfilt_max;j++){
       FitFluxmin[i][j]=99999.;
       DataFluxmin[i][j]=99999.;
     }
     DataFlux_min[i]=99999.;
     FitFlux_min[i]=99999.;
 }

 const Float_t Error_min=1.e-5;
 Float_t SNRmax_data[Nsnrmax]={0.}, snrvalue;
 Int_t found, npass, ntitlelines;
 Int_t DataFlag[Ntot_max][Nfilt_max];
 Int_t Nobs[Nplotmax][Nfilt_max]={0};
 Int_t Ndata[Nplotmax][Nfilt_max]={0},Nfit[Nplotmax][Nfilt_max]={0};
 Int_t Ngraphs[Nplotmax]={0};
 char CID[Nplotmax][nchar], TREEID[Nplotmax][nchar]; 
 char displaytext0[nchar]="not_found"; char displaytext1[nchar]="not_found"; 
 char displaytext2[nchar]="not_found"; char displaytext3[nchar]="not_found";
 char title[Nplotmax][nchar]={""};
 TGraphErrors *flux_vs_epoch[Nplotmax][Nfilt_max]; 
 TGraph *bestfit_vs_epoch_noerrors[Nplotmax][Nfilt_max];
 TGraphErrors *bestfit_vs_epoch[Nplotmax][Nfilt_max];

 Int_t ntree, nsne;

n=-1;   //Initialize plot counter to -1 since arrays start with index 0 
//Loop over trees
for (ntree=0;ntree<nlcpak;ntree++) {
 TTree *snlcpak = (TTree*)f->Get(snlcpak_tree[ntree]);
 if(pchr=strstr(printoptions,"print3")!= NULL) snlcpak->Print();

 if(snlcpak->GetBranchStatus("CCID")==1) snlcpak->SetBranchAddress("CCID",Ccid);
 if(snlcpak->GetBranchStatus("NTOT")==1) snlcpak->SetBranchAddress("NTOT",&Ntot);
 if(snlcpak->GetBranchStatus("IFILTOBS")==1) snlcpak->SetBranchAddress("IFILTOBS",Ifiltobs);
 if(snlcpak->GetBranchStatus("CFILTOBS")==1) snlcpak->SetBranchAddress("CFILTOBS",Cfiltobs);
 if(snlcpak->GetBranchStatus("FLUXCAL")==1) snlcpak->SetBranchAddress("FLUXCAL",Fluxcal_data);
 if(snlcpak->GetBranchStatus("FLUXCAL_ERR")==1) snlcpak->SetBranchAddress("FLUXCAL_ERR",Fluxcal_data_error);
 if(snlcpak->GetBranchStatus("TOBS")==1) snlcpak->SetBranchAddress("TOBS",Tobs_data);
 if(snlcpak->GetBranchStatus("FITCHI2")==1) snlcpak->SetBranchAddress("FITCHI2",Fitchi2);
 if(snlcpak->GetBranchStatus("IFLAGDATA")==1) snlcpak->SetBranchAddress("IFLAGDATA",Iflagdata);
 if(snlcpak->GetBranchStatus("REJECT")==1) snlcpak->SetBranchAddress("REJECT",Reject);
 if(snlcpak->GetBranchStatus("DISPLAYTEXT0")==1) snlcpak->SetBranchAddress("DISPLAYTEXT0",displaytext0);
 if(snlcpak->GetBranchStatus("DISPLAYTEXT1")==1) snlcpak->SetBranchAddress("DISPLAYTEXT1",displaytext1);
 if(snlcpak->GetBranchStatus("DISPLAYTEXT2")==1) snlcpak->SetBranchAddress("DISPLAYTEXT2",displaytext2);
 if(snlcpak->GetBranchStatus("DISPLAYTEXT3")==1) snlcpak->SetBranchAddress("DISPLAYTEXT3",displaytext3);
 if(snlcpak->GetBranchStatus("BAND_PEAKMJD")==1) snlcpak->SetBranchAddress("BAND_PEAKMJD",Band_peakmjd);
 
 Int_t nentries = Int_t(snlcpak->GetEntriesFast());
 cout<<" Number of "<<snlcpak_tree[ntree]<<" Entries (SNe) "<<nentries<<"\n"<<endl;
 Int_t nbytesSK=0, nbSK=0;

 for (Int_t jentry=0; jentry<nentries;jentry++) {
   Int_t ientry = snlcpak->LoadTree(jentry); if (ientry < 0) break;
   nbSK = snlcpak->GetEntry(jentry);   nbytesSK += nbSK;
   if(pchr=strstr(printoptions,"print3")!= NULL) {
      cout<<"  Loading entry "<<jentry<<" with cid = "<<Ccid<<";  Number of observations = "<<Ntot<<endl;
      cout<<"  Header Information"<<endl;
      cout<<"  :"<<displaytext0<<endl;
      cout<<"  :"<<displaytext1<<endl;
      cout<<"  :"<<displaytext2<<endl;
      cout<<"  :"<<displaytext3<<endl;
      cout<<" "<<endl;
   }
   if(Ntot > Ntot_max) {
     cout<<"Warning: Maximum number of observations exceeded: Ntot = "<<Ntot<<" Skipping SN #"<<Ccid<<endl;
   }
   else {
    if( ( SNid==0 && jentry==nentries-1) || (SNid <0 && jentry<NSNe) ||
        ( pchr=strstr(Ccid,SNcid)!= NULL && strlen(Ccid)==strlen(SNcid)) ) {

      if(pchr=strstr(options,"SNR_")!= NULL) {  //SNR cut requested?
	for (k=0;k<Nsnr;k++) {
	  SNRmax_data[k]=0;
        }
        for (j=0;j<Ntot;j++) {  //Find snrmax values
	  if( Iflagdata[j]==1  && Ifiltobs[j]>=Ifilt_min && Ifiltobs[j]<=Ifilt_max) {
	    if(Fluxcal_data_error[j]>Error_min) {
              snrvalue = Fluxcal_data[j]/Fluxcal_data_error[j];
              found=0;
              for (k=0;k<Nsnr;k++) {
                if(snrvalue > SNRmax_data[k] && found==0) {
                  for(l=Nsnr-1;l>k;l--){
		    SNRmax_data[l]=SNRmax_data[l-1];  //move entries down 
                  }
                  SNRmax_data[k]=snrvalue;
                  found=1;
                }
              }
              if(pchr=strstr(printoptions,"debug")!= NULL){
               cout<<j<<","<<Iflagdata[j]<<","<<Fluxcal_data[j]/Fluxcal_data_error[j]<<endl;
                 for (l=0;l<Nsnr;l++) {
                    cout<<"SNRmax #"<<l+1<<" in data = "<<SNRmax_data[l]<<endl;
	         }
	      }
	    }
	  }
	} //snrmax values
        npass=0;
        for (k=0;k<Nsnr;k++) {
         cout<<"SNRmax #"<<k+1<<" in data = "<<SNRmax_data[k]<<endl;
         if(SNRmax_data[k]>=SNRmax[k]) npass++;
        }
        if (npass == Nsnr) {
           cout<<"Candidate passes "<<snroption<<" cuts\n"<<endl;
	}
        else {
           cout<<"Candidate fails "<<snroption<<" cuts\n"<<endl;
        }
      }

      if(npass == Nsnr) {
       n=n+1;     //Increment plot counter
       cout<<" Found matching SN cid "<<Ccid<<" for request "<<SNcid<<": Filling plot #"<<n+1<<"\n"<<endl;
       sprintf(CID[n],"%s",Ccid);
       sprintf(TREEID[n],"%s",snlcpak_tree[ntree]);
       cout<<" Band Peak MJD's are:"<<endl;
       for (j=0;j<Ifilt_max-Ifilt_min+1;j++){
         cout<<" Filter:"<<j<<" "<<Band_peakmjd[j]<<endl;
       }
       //Save graph title
       //ntitlelines=0;
       //sprintf(title[n],"#splitline");
       sprintf(title[n],"   ");
       if(strlen(displaytext0)>0) {
         if (pchr1=strstr(displaytext0,"NULL")== NULL) {
	   //sprintf(title[n],"%s{%s}",title[n],displaytext0);
	   sprintf(title[n],"%s %s",title[n],displaytext0);
	   //ntitlelines++;
	 }
       }
       if (strlen(displaytext1) > 0) {
         if (pchr1=strstr(displaytext1,"NULL")== NULL) {
	   //sprintf(title[n],"%s{%s}",title[n],displaytext1);
	   sprintf(title[n],"%s %s",title[n],displaytext1);
	   //cout<<title[n]<<endl;
	   //ntitlelines++;
	 }
       }
       if (strlen(displaytext2) > 0) {
         if (pchr1=strstr(displaytext2,"NULL")== NULL) {
	   //sprintf(title[n],"%s{%s}",title[n],displaytext2);
	   sprintf(title[n],"%s %s",title[n],displaytext2);
	   //cout<<title[n]<<endl;
	   //ntitlelines++;
	 }
       } 
       if (strlen(displaytext3) > 0) {
         if (pchr1=strstr(displaytext3,"NULL")== NULL) {
	   //sprintf(title[n],"%s{%s}",title[n],displaytext3);
	   sprintf(title[n],"%s %s",title[n],displaytext3);
	   //cout<<title[n]<<endl;
	   //ntitlelines++;
	 }
       } 
       if(pchr=strstr(printoptions,"print1")!= NULL){
	 // cout<<"Saved graph title "<<title[n]<<" with "<<ntitlelines<<" lines"<<"\n"<<endl;
	 cout<<"Saved graph title "<<title[n]<<"\n"<<endl;
       }

       if(pchr=strstr(printoptions,"print2")!= NULL) cout<<"obs: Iflagdata  Reject Tobs  Ifiltobs Cfiltobs Fluxcal_data Fluxcal_error"<<endl;
       for (j=0;j<Ntot;j++){
        if(Ifiltobs[j]<Ifilt_min || Ifiltobs[j]>Ifilt_max) {
          cout<<"Warning! Filter index = "<<Ifiltobs[j]<<" out of allowed range. Skipping entry "<<j<<endl;
        }
        else{
          //Organize data by filternumber
          if (Fluxcal_data[j] > fluxlimit) {
              DataFlag[ Nobs[n][Ifiltobs[j]] ][Ifiltobs[j]] = Iflagdata[j];
              Epoch[ Nobs[n][Ifiltobs[j]] ][Ifiltobs[j]] = Tobs_data[j];
              Flux[ Nobs[n][Ifiltobs[j]] ][Ifiltobs[j]]  = Fluxcal_data[j];
              Flux_error[ Nobs[n][Ifiltobs[j]] ][Ifiltobs[j]]  = Fluxcal_data_error[j];
              Nobs[n][Ifiltobs[j]]++;
              if(Iflagdata[j]==1) { //real data
                 Ndata[n][Ifiltobs[j]]++;
                 DataEpochmin = TMath::Min(DataEpochmin,Tobs_data[j]);
                 DataEpochmax = TMath::Max(DataEpochmax,Tobs_data[j]);
                 DataFluxmax[n][Ifiltobs[j]] = TMath::Max(DataFluxmax[n][Ifiltobs[j]],Fluxcal_data[j]);  //Max in each filter
                 DataFlux_max[n] = TMath::Max(DataFlux_max[n],Fluxcal_data[j]);  //Max for all filters
                 DataFluxMAX = TMath::Max(DataFluxMAX,Fluxcal_data[j]);  //Max for all SNe
                 DataFluxmin[n][Ifiltobs[j]] = TMath::Min(DataFluxmin[n][Ifiltobs[j]],Fluxcal_data[j]);  //Min in each filter
                 DataFlux_min[n] = TMath::Min(DataFlux_min[n],Fluxcal_data[j]);  //Min for all filters
                 DataFluxMIN = TMath::Min(DataFluxMIN,Fluxcal_data[j]);  //Min for all SNe
                 if(pchr=strstr(options,"debug")!= NULL) cout<<DataFluxmin[n][Ifiltobs[j]]<<","<<DataFlux_min[n]<<","<<DataFluxMIN<<endl;
	      }
              else if(Iflagdata[j]==0) { //fit data
                 Nfit[n][Ifiltobs[j]]++;
                 FitEpochmin = TMath::Min(FitEpochmin,Tobs_data[j]);
                 FitEpochmax = TMath::Max(FitEpochmax,Tobs_data[j]);
                 FitFluxmax[n][Ifiltobs[j]] = TMath::Max(FitFluxmax[n][Ifiltobs[j]],Fluxcal_data[j]);  //Max in each filter
                 FitFlux_max[n] = TMath::Max(FitFlux_max[n],Fluxcal_data[j]);  //Max for all filters
                 FitFluxMAX = TMath::Max(FitFluxMAX,Fluxcal_data[j]);  //Max for all SNe
                 FitFluxmin[n][Ifiltobs[j]] = TMath::Min(FitFluxmin[n][Ifiltobs[j]],Fluxcal_data[j]);  //Min in each filter
                 FitFlux_min[n] = TMath::Min(FitFlux_min[n],Fluxcal_data[j]);  //Min for all filters
                 FitFluxMIN = TMath::Min(FitFluxMIN,Fluxcal_data[j]);  //Min for all SNe
                 if(pchr=strstr(options,"debug")!= NULL) cout<<FitFluxmin[n][Ifiltobs[j]]<<","<<FitFlux_min[n]<<","<<FitFluxMIN<<endl;
              }
             if(pchr=strstr(printoptions,"print2")!= NULL) 
                cout<<j<<":    "<<Iflagdata[j]<<",        "<<Reject[j]<<",  "<<Tobs_data[j]<<",     "<<Ifiltobs[j]<<",    "<<Cfiltobs[j]<<",     "<<Fluxcal_data[j]<<",     "<<Fluxcal_data_error[j]<<endl;
          }
	} //endif on Ifilt check
       } //endfor on observations

       Int_t ifilter; 
       //Calculate Fluxmax and min arrays depending on plot options
       if(pchr=strstr(options,"data")!= NULL) {  //data only
	for (i=0;i<Nfilters;i++) {
	  ifilter = IfilterId[i];
	  Fluxmax[n][ifilter] = DataFluxmax[n][ifilter];
	  Fluxmin[n][ifilter] = DataFluxmin[n][ifilter];
	} 
	Epochmin = DataEpochmin;
	Epochmax = DataEpochmax;
	Flux_max[n] = DataFlux_max[n];
        FluxMAX = DataFluxMAX;
	Flux_min[n] = DataFlux_min[n];
        FluxMIN = DataFluxMIN;
       }
       else if(pchr=strstr(options,"fit")!= NULL) {  //fit only
	for (i=0;i<Nfilters;i++) {
	  ifilter = IfilterId[i];
	  Fluxmax[n][ifilter] = FitFluxmax[n][ifilter];
          Fluxmin[n][ifilter] = FitFluxmin[n][ifilter];
	} 
	Flux_max[n] = FitFlux_max[n];
	Epochmin = FitEpochmin;
	Epochmax = FitEpochmax;
        FluxMAX = FitFluxMAX;
        Flux_min[n] = FitFlux_min[n];
        FluxMIN = FitFluxMIN;
       }
       else {
	for (i=0;i<Nfilters;i++) {
	  ifilter = IfilterId[i];
	  Fluxmax[n][ifilter] = TMath::Max(DataFluxmax[n][ifilter],FitFluxmax[n][ifilter]);
	  Fluxmin[n][ifilter] = TMath::Min(DataFluxmin[n][ifilter],FitFluxmin[n][ifilter]);
	} 
	Flux_max[n] = TMath::Max(DataFlux_max[n],FitFlux_max[n]);
	Flux_min[n] = TMath::Min(DataFlux_min[n],FitFlux_min[n]);
	Epochmin = DataEpochmin; //use range in data
	Epochmax = DataEpochmax;
        FluxMAX = TMath::Max(DataFluxMAX,FitFluxMAX);
        FluxMIN = TMath::Min(DataFluxMIN,FitFluxMIN);
       }
       if(pchr=strstr(printoptions,"print2")!= NULL) {
        for (i=0;i<Nfilters;i++) {
	  cout<<"Maximum flux for filter "<<IfilterId[i]<<" = "<<Fluxmax[n][IfilterId[i]]<<endl;
	  cout<<"Minimum flux for filter "<<IfilterId[i]<<" = "<<Fluxmin[n][IfilterId[i]]<<endl;
	}
	cout<<"Maximum flux over filters"<<" = "<<Flux_max[n]<<endl;
	cout<<"Minimum flux over filters"<<" = "<<Flux_min[n]<<endl;
        cout<<"Epochmin = "<<Epochmin<<"; Epochmax = "<<Epochmax<<"\n"<<endl;
       }
       //Copy arrays into temporary vectors and fill TGraphs
       Float_t Epoch_temp[Ntot_max], Flux_temp[Ntot_max], Fluxerr_temp[Ntot_max],Epocherr_temp[Ntot_max]={0};
       Float_t FitEpoch_temp[Ntot_max], FitFlux_temp[Ntot_max], FitFluxerr_temp[Ntot_max],FitEpocherr_temp[Ntot_max]={0};
       Int_t ndata,nfit;

       for (i=0;i<Nfilters;i++) {
	ifilter = IfilterId[i];
	if(pchr=strstr(printoptions,"print3")!= NULL) 
          cout<<"Check of Array assignment by filter: "<<ifilter<<" Nobs = "<<Nobs[n][ifilter]<<endl;
	if(Nobs[n][ifilter]>0){
	  Ngraphs[n]=Ngraphs[n]+1;
	  ndata=-1; nfit=-1;
	  for ( j=0;j<Nobs[n][ifilter];j++){
	    if(DataFlag[j][ifilter]==1) {
	      ndata=ndata+1;
	      Epoch_temp[ndata] = Epoch[j][ifilter];
              if (pchr=strstr(yscale,"mag")!= NULL) {
	        Flux_temp[ndata] = 2.5*log10(Flux[j][ifilter])-zeropt;
	        Fluxerr_temp[ndata] = 2.5/log(10)*Flux_error[j][ifilter]/Flux[j][ifilter];
	        if(pchr=strstr(printoptions,"print3")!= NULL) 
                 cout<<"ifilter="<<ifilter<<", Epoch="<<Epoch_temp[ndata]<<", Data Magnitude="<<-Flux_temp[ndata]<<endl;
              }
              else {
	        Flux_temp[ndata] = Flux[j][ifilter];
	        Fluxerr_temp[ndata] = Flux_error[j][ifilter];
	        if(pchr=strstr(printoptions,"print3")!= NULL) 
                 cout<<"ifilter="<<ifilter<<", Epoch="<<Epoch_temp[ndata]<<", Data Flux="<<Flux_temp[ndata]<<endl;
              }
	    }
	    else if(DataFlag[j][ifilter]==0) {
	      nfit=nfit+1;
	      FitEpoch_temp[nfit] = Epoch[j][ifilter];
              if (pchr=strstr(yscale,"mag")!= NULL) {
	        FitFlux_temp[nfit] = 2.5*log10(Flux[j][ifilter])-zeropt;
	        FitFluxerr_temp[nfit] = 2.5/log(10)*Flux_error[j][ifilter]/Flux[j][ifilter];
	        if(pchr=strstr(printoptions,"print3")!= NULL) 
		 cout<<"ifilter="<<ifilter<<", Fit Epoch="<<FitEpoch_temp[nfit]<<", Fit Magnitude="<<-FitFlux_temp[nfit]<<endl; 
              }
              else {
	        FitFlux_temp[nfit] = Flux[j][ifilter];
	        FitFluxerr_temp[nfit] = Flux_error[j][ifilter];
	        if(pchr=strstr(printoptions,"print3")!= NULL) 
                 cout<<"ifilter="<<ifilter<<", Fit Epoch="<<FitEpoch_temp[nfit]<<", Fit Flux="<<FitFlux_temp[nfit]<<endl;
              }
	    }
	  }
	  if(ndata+1!=Ndata[n][ifilter]) cout<<"Ooops, warning: mismatch in number of data points counted "<<ndata<<","<<Ndata[n][ifilter]<<endl;
	  flux_vs_epoch[n][i] = new TGraphErrors(Ndata[n][ifilter],Epoch_temp,Flux_temp,Epocherr_temp,Fluxerr_temp);
	  if(nfit+1!=Nfit[n][ifilter]) cout<<"Ooops, warning: mismatch in number of data points counted "<<nfit<<","<<Nfit[n][ifilter]<<endl;
	  bestfit_vs_epoch[n][i] = new TGraphErrors(Nfit[n][ifilter],FitEpoch_temp,FitFlux_temp,FitEpocherr_temp,FitFluxerr_temp);
	  bestfit_vs_epoch_noerrors[n][i] = new TGraph(Nfit[n][ifilter],FitEpoch_temp,FitFlux_temp);
          if(pchr=strstr(printoptions,"print2")!= NULL) 
            cout<<"Saved "<<Ndata[n][ifilter]<<" data points and "<<Nfit[n][ifilter]<<" fit points for filter #"<<ifilter<<" in plot #"<<n+1<<endl;
	  /*
	  if (Ndata[n][ifilter]==0 && pchr=strstr(options,"fit")== NULL)  {
	    cout<<"Note: No data points found for SN #"<<CID[n]<<" in filter "<<ifilter<<endl;
	  }
	  if (Nfit[n][ifilter]==0 && pchr=strstr(options,"data")== NULL)  {
            cout<<"Note: No fit points found for SN #"<<CID[n]<<" in filter "<<ifilter<<endl;
	  }
          */
	}//if-Nobs>0
       }
       cout<<"Saved "<<Ngraphs[n]<<" (non-empty) graphs for SN #"<<CID[n]<<" in tree "<<ntree<<" ("<<TREEID[n]<<")\n"<<endl;

      } //endif on npass match
    } //endif on SN match
   } //Check on Ntot       
 } //loop over entries 
}//end of loop over trees

 Int_t Nplot, NSN_tree;
 Nplot=n+1;
 cout<<"Total number of graphs saved = "<<Nplot<<"\n"<<endl;
 NSN_tree=Nplot/nlcpak;

 if(nlcpak>1) {  //check that Nplot is a multiple of nlcpak
   if(Nplot%nlcpak!=0) cout<<"Warning: mismatch ?. Number of graphs not equal to (number of trees) x (number selected SNe per tree)\n"<<endl;
 }

//Check for FITRES tree
// Key code
 char fitresCcid[nchar];
 Int_t nfound, zfound;
 char z_type[nchar]; char z_photo[nchar]="photo-z"; char z_spec[nchar]="spec-z";
 Float_t z, zerror, z_value[Nplotmax]={-1.},dz_value[Nplotmax]={0.};
 
 if(pchr=strstr(printoptions,"debug")!= NULL) f->GetListOfKeys()->Print();
 TIter next(f->GetListOfKeys());
 TKey *key; TTree *tree;
 while ((key=(TKey*)next())) {
   if(strstr("TTree",key->GetClassName())!= NULL) {
     if(strstr("FITRES",key->GetName())!= NULL) {
       cout<<"FITRES tree available: adding z information "<<endl;
       fitres=(TTree*)f->Get(key->GetName());
       if(pchr=strstr(printoptions,"print3")!= NULL) fitres->Print();
       if(fitres->GetBranchStatus("CCID")==1) fitres->SetBranchAddress("CCID",fitresCcid);
       zfound=0;
       if(fitres->GetBranchStatus("ZPHOT")==1) {
          fitres->SetBranchAddress("ZPHOT",&z);
          if(fitres->GetBranchStatus("ZPHOTERR")==1) {
             fitres->SetBranchAddress("ZPHOTERR",&zerror);
             zfound=1;
          }
          sprintf(z_type,"%s",z_photo);
       }
       else if(fitres->GetBranchStatus("z")==1) {
          fitres->SetBranchAddress("z",&z);
          if(fitres->GetBranchStatus("zERR")==1) {
             fitres->SetBranchAddress("zERR",&zerror);
             zfound=1;
          }
          sprintf(z_type,"%s",z_spec);
       }
       else {
	 cout<<"Branch address for z not found"<<endl;
       }
       if(fitres->GetBranchStatus("CCID")==1 && zfound==1) {
         nentries = Int_t(fitres->GetEntriesFast());
         cout<<" Number of "<<key->GetName()<<" Entries (SNe) "<<nentries<<"\n"<<endl;
         Long64_t nbytes = 0, nb = 0;
         jentry=0;
         nfound=0;
         while (jentry<nentries && nfound<Nplot) {
           ientry = fitres->LoadTree(jentry);
           if (ientry < 0) break;
           nb = fitres->GetEntry(jentry);   nbytes += nb;
           if(pchr1=strstr(printoptions,"print3")!= NULL) cout<<"CID="<<fitresCcid<<", "<<z_type<<"= "<<z<<" +/- "<<zerror<<endl;
           n=0;
           found=0;          
	   while (n<Nplot && found==0) {  //search for matching CID
             if(pchr1=strstr(CID[n],fitresCcid)!= NULL && strlen(CID[n])==strlen(fitresCcid)) {
	       if(pchr2=strstr(printoptions,"print1")!= NULL) cout<<"CID match! Assigning z value of "<<z<<" to SN "<<CID[n]<<endl;
               z_value[n]=z;
	       dz_value[n]=zerror;
               found=1;
             } //cid match 
             n++;
	   } //cid search
           if (found==1) nfound++;
           jentry++;
         }//entries
       }
       else {
         cout<<" Ntuple variables CCID and/or Z and/or ZERR are not available\n"<<endl;
       }
     } //fitres check
   } //loop over trees
 } //loop over keys
 f->Close();


//Draw canvases
 Float_t ymax; Float_t ymin; Float_t ysize;
 Float_t ymaxa,ymaxb,ymaxc; Float_t ymina,yminb,yminc;
 Int_t imax,imin; Float_t xmin, xmax;
 Float_t big_margin=0.1;  Float_t xbig_margin=0.5; Float_t small_margin=0.05;  Float_t med_margin=0.25;
 Float_t vert_scale_factor=1.10;
 Float_t x_lo=0.2; Float_t y_lo=0.15; Float_t x_hi=0.9; Float_t y_hi=0.99;
 Float_t y_4=0.28; Float_t y_3=0.18; Float_t y_2=0.09; Float_t y_1=0.05; Float_t y_fit=0.23;
 Float_t x_vlong=0.7; Float_t x_long=0.65; Float_t x_med=0.45; Float_t x_short=0.4, x_vshort=0.3;
 Float_t title_xlo=0.02, title_xhi=0.98, title_ylo=0.96, title_yhi=0.999, pad_lo=.01, pad_xhi=0.99;
 Float_t ytitleoffset=1.17;
 Float_t lgndtextscale=0.8; //(scale factor for textsize in legends)
 Float_t labeltextscale=0.8; //(scale factor for textsize on axes)
 Int_t lgndtextfont=textfont;
 const Int_t nlgnd=100;
 char lgnd_title[nlgnd]="";  
 char xtitle[nchar],ytitle[nchar];
 Float_t screen_x,screen_y;

 char extraplotopts[nchar];
 if ((pchr=strstr(scaleoption,"multi")) != NULL) {
    sprintf(extraplotopts,"_multiscale");
 }
 else if ((pchr=strstr(scaleoption,"one")) != NULL) {
    sprintf(extraplotopts,"_onescale");
 }
 else {
   sprintf(extraplotopts,"");
 }
 if((pchr=strstr(layout,"column")) != NULL) {
   sprintf(extraplotopts,"%s_cformat",extraplotopts);
 }

 sprintf(xtitle,"T-T_{0} (days)");
 if (pchr=strstr(yscale,"mag")!= NULL) {
   sprintf(ytitle,"Magnitude");
   sprintf(extraplotopts,"%s_mag",extraplotopts);
 }
 else {
   sprintf(ytitle,"Flux");
 }

 if(pchr=strstr(outputoption,"one")!= NULL) {
   if(pchr1=strstr(outputoption,"snana")!= NULL) {
     sprintf(outputfile,"%s.%s",infile,filetype);
   }
   else {
     sprintf(outputfile,"%s%s.%s",infile,extraplotopts,filetype);
   }
   sprintf(openoutputfile,"%s[",outputfile);
   sprintf(closeoutputfile,"%s]",outputfile);
 }
// Setup scale limits if requested
 Float_t fluxlim_lo, fluxlim_hi, maglim_lo, maglim_hi;
 if(strstr(ylimits,"auto")== NULL) {       //fixed scales requested
     if(strstr(yscale,"flux")!= NULL) {  //flux
       fluxlim_lo=ylim_lo;
       fluxlim_hi=ylim_hi;
     }
     else {                                   //mag (will be -ve, so flip values)
       maglim_lo=ylim_hi;
       maglim_hi=ylim_lo;
     }
 }
 Int_t const  Ncanvasmax = Nplotmax; 
 Int_t Nrows,Ncolmax=4;
 Float_t screensize=1000, row_aspectratio=2, col_aspectratio=3;
 Int_t plotid; Int_t Ncolumns;
 TCanvas *c[Ncanvasmax]; 
 char canvas[Ncanvasmax][2]; char canvas_title[Ncanvasmax][nchar]; 
 Int_t isn,itree;

 //Loop over trees and make canvases
 //for (n=0;n<Nplot;n++) { //old loop over plots
 for (isn=0;isn<NSN_tree;isn++) {  //loop first over SNe to keep plots of single SN together
  for (itree=0;itree<nlcpak;itree++) { //next over trees
   n=itree*NSN_tree+isn;                       //compute plot number
   //cout<<isn<<","<<itree<<","<<n<<endl;
   sprintf(canvas[n],"c%d",n+1);
   sprintf(canvas_title[n],"LC_%s",TREEID[n]);
   cout<<"Drawing canvas with "<<Ngraphs[n]<<" graphs for tree "<<TREEID[n]<<endl;

   if(Ngraphs[n]%Ncolmax==0){
    Nrows = Ngraphs[n]/Ncolmax;
   }
   else{
    Nrows = Ngraphs[n]/Ncolmax+1;
   }
   Ncolumns=(long)((double)Ngraphs[n]/(double)Nrows+0.5);
   if(pchr=strstr(printoptions,"debug")!= NULL) cout<<"Nrows,Ncolumns = "<<Nrows<<","<<Ncolumns<<endl;

   if((pchr=strstr(layout,"row")) != NULL) {
    screen_x=screensize;
    screen_y=screen_x/row_aspectratio;
   }
   else {
    screen_y=screensize;   //column layout, flip rows and columns
    screen_x=screen_y/col_aspectratio;
    Int_t isave=Nrows;
    Nrows = Ncolumns;
    Ncolumns = isave;
   }

   Int_t padsep=0.0; 

   char plottype[nchar];
   sprintf(plottype,"");
   if(pchr=strstr(options,"data")!= NULL) sprintf(plottype,"_data");
   if(pchr=strstr(options,"fit")!= NULL) sprintf(plottype,"_fit");

   c[n] = new TCanvas(canvas[n],canvas_title[n],0,0,screen_x,screen_y);  //laptop-screen
   //Set global title
   TPaveLabel* globaltitle = new TPaveLabel(title_xlo,title_ylo,title_xhi,title_yhi,title[n]);
   globaltitle->Draw();

   TPad* graphPad = new TPad("Graphs","Graphs",pad_lo,pad_lo,pad_xhi,title_ylo);
   graphPad->Draw();
   graphPad->cd();
   if((pchr=strstr(options,"merge")) != NULL) {
    x_hi=0.999; y_hi=0.999;
    //c[n]->Divide(Ncolumns,Nrows,padsep,padsep);
    graphPad->Divide(Ncolumns,Nrows,padsep,padsep);
   }
   else {
    //c[n]->Divide(Ncolumns,Nrows);
    graphPad->Divide(Ncolumns,Nrows);
   }

   plotid=0;
   for (i=0;i<Ngraphs[n];i++) {
     plotid++;
     //c[n]->cd(plotid)->SetGrid(1);
     graphPad->cd(plotid)->SetGrid(1);
   }

   plotid=0;
   for (i=0;i<Nfilters;i++) {
    if(Nobs[n][IfilterId[i]]>0){
     plotid++;
     if ((pchr=strstr(scaleoption,"multi")) != NULL) {  //select which max,min
        ymax=Fluxmax[n][IfilterId[i]];
        ymin=Fluxmin[n][IfilterId[i]];
     }
     else if ((pchr=strstr(scaleoption,"one")) != NULL) {
        ymax=FluxMAX; 
        ymin=FluxMIN; 
     }
     else {
        ymax=Flux_max[n]; 
        ymin=Flux_min[n]; 
     }
     if(pchr=strstr(options,"debug")!= NULL)  cout<<"ymax, ymin in Flux scale ="<<ymax<<","<<ymin<<endl;
     if (pchr=strstr(yscale,"mag")!= NULL) {     //convert to -ve magnitudes
        if(pchr1=strstr(ylimits,"max")!= NULL) {
            ymax = - maglim_hi;
        }
        else {        
            ymax=2.5*log10(ymax)-zeropt;
        }
        if(pchr1=strstr(ylimits,"min")!= NULL) {
            ymin = - maglim_lo;
        }
        else {
            ymin=2.5*log10(ymin)-zeropt;
        }
     }
     else {        //flux option
       if(pchr1=strstr(options,"posflux")!= NULL) ymin=0.;  //force scale to zero for posflux option
       if(pchr1=strstr(ylimits,"max")!= NULL) ymax = fluxlim_hi;  //force scale
       if(pchr1=strstr(ylimits,"min")!= NULL) ymin = fluxlim_lo;  //force scale
     }
     ysize=ymax-ymin;
     if(pchr1=strstr(ylimits,"auto")!= NULL) {  //add margins for any autoscale
       ymax=ymax+big_margin*(ysize); ymin=ymin-small_margin*(ysize);
     }
     elseif(pchr1=strstr(ylimits,"minonly")!= NULL) {
       ymax=ymax+big_margin*(ysize);
     }
     elseif(pchr1=strstr(ylimits,"maxonly")!= NULL) {
       ymin=ymin-small_margin*(ysize);
     }
     if(pchr=strstr(printoptions,"print1")!= NULL) {
       cout<<ytitle<<" scale for plot of filter #"<<IfilterId[i]<<"="<<ymin<<" to "<<ymax<<endl;
     }   
     xmax=Epochmax+small_margin*(Epochmax-Epochmin); xmin=Epochmin-small_margin*(Epochmax-Epochmin);

     //c[n]->cd(plotid);
     graphPad->cd(plotid);
     sprintf(lgnd_title," %s band  SN %s   ",FilterId[i],CID[n]);
     if(z_value[n] > 0) {
         sprintf(lgnd_title,"%s z=%.2f#pm%.2f",lgnd_title,z_value[n],dz_value[n]);
     }
     if(pchr=strstr(printoptions,"print2")!= NULL) cout<<"Drawing "<<lgnd_title<<" graph"<<endl;
     TLegend* lega = new TLegend(x_hi-x_long,y_hi-y_2,x_hi,y_hi);
     lega->SetNColumns(2);     
     lega->SetTextFont(lgndtextfont);
     lega->SetHeader(lgnd_title);
     lega->SetTextSize(lgndtextscale*textsize); 
     if(pchr=strstr(options,"fit")== NULL && Ndata[n][IfilterId[i]] > 0) {  //data requested and data points found
       flux_vs_epoch[n][i]->GetHistogram()->SetMaximum(ymax);
       flux_vs_epoch[n][i]->GetHistogram()->SetMinimum(ymin);
       flux_vs_epoch[n][i]->GetXaxis()->SetLimits(xmin,xmax);
       flux_vs_epoch[n][i]->SetMarkerColor(kBlack);
       flux_vs_epoch[n][i]->SetTitle("");
       flux_vs_epoch[n][i]->GetXaxis()->SetTitleFont(textfont); flux_vs_epoch[n][i]->GetXaxis()->SetTitleSize(textsize);
       flux_vs_epoch[n][i]->GetXaxis()->SetLabelFont(textfont); flux_vs_epoch[n][i]->GetXaxis()->SetLabelSize(labeltextscale*textsize);
       flux_vs_epoch[n][i]->GetXaxis()->SetTitle(xtitle); flux_vs_epoch[n][i]->GetXaxis()->CenterTitle();
       flux_vs_epoch[n][i]->GetYaxis()->SetTitleFont(textfont); flux_vs_epoch[n][i]->GetYaxis()->SetTitleSize(textsize);
       flux_vs_epoch[n][i]->GetYaxis()->SetLabelFont(textfont); flux_vs_epoch[n][i]->GetYaxis()->SetLabelSize(labeltextscale*textsize);
       //if(plotid%Ncolumns==1) {
         flux_vs_epoch[n][i]->GetYaxis()->SetTitle(ytitle); 
         flux_vs_epoch[n][i]->GetYaxis()->SetTitleOffset(ytitleoffset);flux_vs_epoch[n][i]->GetYaxis()->CenterTitle(); 
       //}
       //if(pchr2=strstr(options,"data")== NULL  &&  Nfit[n][IfilterId[i]] > 0) { //fit already plotted?
       //  flux_vs_epoch[n][i]->Draw("p e1");
       //} 
       //else {
         flux_vs_epoch[n][i]->Draw("a p e1");
       //}
       lega->AddEntry(flux_vs_epoch[n][i],"Data","P");
     }
     if(pchr=strstr(options,"redline")!= NULL) {
       //cd(plotid)->Update();
       //TLine *redline=new TLine(cd(plotid)->GetUxmin(),0.0,cd(plotid)->GetUxmax(),0.0);
       if(pchr1=strstr(options,"print3")!= NULL) {
	 cout<<"Adding red line at y=0 for "<<xmin<<" < x < "<<xmax<<"\n"<<endl;
       }
       TLine *redline=new TLine(xmin,0.0,xmax,0.0);
       redline->SetLineColor(kRed);
       redline->Draw();
     }
     if(pchr=strstr(options,"data")== NULL &&  Nfit[n][IfilterId[i]] > 0) {  //fit requested and fit points found
       bestfit_vs_epoch[n][i]->GetHistogram()->SetMaximum(ymax);
       bestfit_vs_epoch[n][i]->GetHistogram()->SetMinimum(ymin);     
       bestfit_vs_epoch_noerrors[n][i]->SetLineColor(kRed);
       bestfit_vs_epoch[n][i]->SetTitle(""); 
       bestfit_vs_epoch[n][i]->GetXaxis()->SetTitleFont(textfont); bestfit_vs_epoch[n][i]->GetXaxis()->SetTitleSize(textsize);
       bestfit_vs_epoch[n][i]->GetXaxis()->SetLabelFont(textfont); bestfit_vs_epoch[n][i]->GetXaxis()->SetLabelSize(labeltextscale*textsize);       
       bestfit_vs_epoch[n][i]->GetXaxis()->SetTitle(xtitle); bestfit_vs_epoch[n][i]->GetXaxis()->CenterTitle();
       bestfit_vs_epoch[n][i]->SetFillColor(44);
       bestfit_vs_epoch[n][i]->GetYaxis()->SetTitleFont(textfont); bestfit_vs_epoch[n][i]->GetYaxis()->SetTitleSize(textsize);
       bestfit_vs_epoch[n][i]->GetYaxis()->SetLabelFont(textfont); bestfit_vs_epoch[n][i]->GetYaxis()->SetLabelSize(labeltextscale*textsize);
       if(pchr2=strstr(options,"fit")== NULL  && Ndata[n][IfilterId[i]] > 0) { //data already plotted ?
         bestfit_vs_epoch[n][i]->Draw("c3");
         bestfit_vs_epoch_noerrors[n][i]->Draw("c");
       } 
       else {
         bestfit_vs_epoch[n][i]->Draw("a3");
         bestfit_vs_epoch_noerrors[n][i]->Draw("c");          
       }
       //if(plotid%Ncolumns==1) {
         bestfit_vs_epoch[n][i]->GetYaxis()->SetTitle(ytitle); 
         bestfit_vs_epoch[n][i]->GetYaxis()->SetTitleOffset(ytitleoffset); bestfit_vs_epoch[n][i]->GetYaxis()->CenterTitle(); 
       //}
       lega->AddEntry(bestfit_vs_epoch[n][i],"Best Fit","F");
       //lega->AddEntry(bestfit_vs_epoch_noerrors[n][i],"Best Fit","L");
     }
     if(pchr=strstr(options,"fit")== NULL && Ndata[n][IfilterId[i]] > 0) {  //data requested and data points found
       flux_vs_epoch[n][i]->Draw("p e1");  //replot data points so they are on top
      }
     //TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),ymin,ymax,510,"+L");
     lega->Draw();
    } //Nobs-check
   }

   if(pchr=strstr(outputoption,"one")!= NULL) {
    if(n==0) c[0]->Print(openoutputfile);  
   }
   else {
    sprintf(outputfile,"%s_%s_CID%s%s%s.%s",infile,TREEID[n],CID[n],plottype,extraplotopts,filetype);
   }
   c[n]->Print(outputfile);

  } //end of loop over trees
 } //end of loop over SNe per tree

 if(pchr=strstr(outputoption,"one")!= NULL) c[0]->Print(closeoutputfile);  

} //end-int
