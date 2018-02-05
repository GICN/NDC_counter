# ---NDC COUNTER------
#
#--Authors: Helene Benveniste (IPSL) & Olivier Boucher (LMD, IPSL)
#--November 2017
#--Probabilistic version
# 
#--to be run in python 2.7
#--execfile('counter_edited_final.py')

print("Launching model")

#
#--import required python packages
import numpy as np
import warnings
from random import random
import matplotlib.pyplot as plt
import math
from math import log, exp
from scipy import interpolate
#
#--key to activate GDP revision according to past values
#revise_GDP=False
revise_GDP=True
#
#--key to activate more verbose during execution
verbose=False
if not verbose: warnings.filterwarnings("ignore")
#
#---------------------READING EMISSIONS DATA--------------------------------------------
#
#--range of future years
future=range(2015,2032)
#
#--opening emission files
path='inputs/'
fid=open(path+'emission_data_EDGAR_FAO.csv','r')
fidLUC=open(path+'emissions_LUC_Julia.csv','r')
fidLUCc=open(path+'list_of_country_codes.csv','r')
fidU=open(path+'emission_data_UNFCCC.csv','r')
fidSSP=open(path+'5regions_SSP_all_CO2eqKyotoGases.csv','r')  
nblineheader=3
nblineheaderU=1
nbyear=5 #--1990-1995-2000-2005-2010 years for which emission data values are provided
#--2012 is not taken into account because EDGAR does not provide fully reliable data for 2012
#--list of greenhouse gases considered
emtypelist=['co2','ch4','n2o','hfc','pfc','sf6','co2lu']  #--in the order of the input file
#
#--All emission values, when not otherwise specified, are expressed in ktCO2eq.
#
#--GWP values organized in a dictionary (gwp)
gwpch4=21. #--SAR value
#gwpch4=25. #--AR4 value
gwpn2o=310.
gwpsf6=23900.
gwp={'co2':1.0, 'ch4':gwpch4,'n2o':gwpn2o,'hfc':1.0,'pfc':1.0,'sf6':gwpsf6,'co2lu':1.0} 
#
#--creating empty dictionaries
emi={}    #--emissions from EDGAR database
emiLUC={} #--emissions LULUCF from BLUE database
emiU={}   #--emissions LULUCF from UNFCCC
emiSSP={} #--emissions projections for SSP on 5 world regions
World=[]  #--list of world countries + transportation sectors
#
#--reading through the EDGAR data-----------------------------------------------------------------
nbline=0
for fline in fid:
  nbline=nbline+1
  if nbline > nblineheader:
    #--reading and splitting line
    data=fline.split(';')
    country=data[0]
    emicountry={}
    #--emission dictionary
    for off,emtype in zip(range(len(emtypelist)),emtypelist):
      emission={}
      for ind in range(0,nbyear):
        year=1990+5*ind
        if data[1+off*nbyear+ind]<>'': 
          emission[year] =float(data[1+off*nbyear+ind])
        else: 
          if verbose: print 'Attention',country, emtype, year, ' missing'
          emission[year] =0.0      
      emicountry[emtype] =emission
    emi[country]=emicountry
    World.append(country)
#
#--reading through the BLUE data for LULUCF---------------------------------------------------------
emiLUCc={}
codeLUC=[]
nbline=0
for fline in fidLUC:
  nbline=nbline+1
  if nbline > nblineheaderU:
    #--reading and splitting line
    data=fline.split(';')
    code=data[0]
    emicode={}
    #--emission dictionary
    for ind in range(0,nbyear):
      year=1990+5*ind
      if data[1+ind]<>'':
        emicode[year] =float(data[1+ind])
      else:
        if verbose: print 'Attention',code, year, ' missing'
        emicode[year] =0.0
    emiLUCc[code]=emicode
    codeLUC.append(code)
#
codelist={}
WorldLUC=[]
nbline=0
for fline in fidLUCc:
  nbline=nbline+1
  if nbline > nblineheaderU+1:
    #--attributing country names to country codes
    data=fline.split(';')
    code=data[1]
    name=data[2].replace("\r\n","")
    codelist[name]=code
    WorldLUC.append(name)
for name in WorldLUC:
  if codelist[name] in codeLUC and codelist[name]<>'0':
    emiLUC[name]=emiLUCc[codelist[name]]
#
#--replacing LUCF emissions from FAO by LUCF emissions from BLUE, and converting from gC to ktCO2.
for country in World:
  for year in range(1990,2015,5):
    if country in WorldLUC and codelist[country]<>'0' and codelist[country] in codeLUC:
      emi[country]['co2lu'][year]=emiLUC[country][year]/1.e9*44./12.
    else:
      emi[country]['co2lu'][year]=0.0
#
del emiLUCc, codeLUC
#
#--fraction of Liu taken for Chinese emissions
#--Liu Z et al. 2015 Reduced carbon emission estimates from fossil fuel combustion and cement 
#--production in China Nature 524 335-338 (doi: 10.1038/nature14677)
#--0= EDGAR, 1= Liu
alphaLiu=0.5
#
#--Dealing with China emissions and Liu correction factors
#--Liu correction factors are read on the graph of Liu's paper
Liu={1990:-0.05,1995:-0.1,2000:-0.06,2005:-0.12,2010:-0.13}
#
#--correct Chinese emissions using mix of EDGAR and Liu
for year in range(1990,2015,5):
  emi['China']['co2'][year]=emi['China']['co2'][year]*(1.+alphaLiu*Liu[year])
#
#
#--reading through the UNFCCC data---------------------------------------------------------------------
nbline=0
for fline in fidU:
  nbline=nbline+1
  if nbline > nblineheaderU:
    #--reading and splitting line
    data=fline.split(';')
    country=data[0]
    emiUcountry={}
    #--UNFCCC emission dictionary
    for ind in range(0,nbyear):
      year=1990+5*ind
      if data[1+ind]<>'':
        emiUcountry[year] =float(data[1+ind])*1000.
      else:
        emiUcountry[year] =0.0
    emiU[country]=emiUcountry
del emiUcountry
#
#--reading through the SSP regions data-----------------------------------------------------------------
#
SSPmodels=['AIM/CGE','IMAGE','WITCH-GLOBIOM','MESSAGE-GLOBIOM','REMIND-MAGPIE']
SSPregions=['ASIA','REF','MAF','LAM','OECD']
SSPscenarios=['SSP1','SSP2','SSP3','SSP4','SSP5']
SSPflavours=['26','34','45','60','Baseline']
#
nbline=0
nbyearSSP=len([2005]+range(2010,2110,10))
for fline in fidSSP:
  nbline=nbline+1
  if nbline > 1:
    #--reading and splitting line
    data=fline.split(';')
    #--finding out the region
    regionfound=False
    for SSPregion in SSPregions: 
      if SSPregion in data[0]: 
        region=SSPregion
        regionfound=True
    #--finding out the model
    modelfound=False
    for SSPmodel in SSPmodels: 
      if SSPmodel in data[1]: 
        model=SSPmodel
        modelfound=True
    #--finding out the scenario
    scenfound=False
    for SSPscen in SSPscenarios: 
      if SSPscen in data[1]: 
        scen=SSPscen
        scenfound=True
    #--finding out the flavour of the scenario
    flavourfound=False
    for SSPflavour in SSPflavours: 
      if SSPflavour in data[1]: 
        flavour=SSPflavour
        flavourfound=True
    if not regionfound or not scenfound or not flavourfound or not modelfound: 
      print 'problem reading SSP emission files', data[1]
      quit()
    if region  not in emiSSP.keys(): emiSSP[region]={}
    if model   not in emiSSP[region].keys(): emiSSP[region][model]={}
    if scen    not in emiSSP[region][model].keys(): emiSSP[region][model][scen]={}
    if flavour not in emiSSP[region][model][scen].keys(): emiSSP[region][model][scen][flavour]={}
    #--emission dictionary
    for ind in range(0,nbyearSSP):
      if ind==0:
        year=2005
      else:
        year=2000+10*ind
      if data[4+ind]<>'':
        emiSSP[region][model][scen][flavour][year]=float(data[4+ind])
      else:
        if verbose: print 'Attention',region, scen, flavour, year, ' missing'
        emiSSP[region][scen][model][flavour][year]=0.0
#
#--------------------------------------------------------------------------------------------------------
#
#--Dealing with country groupings
#
#--EU 28 countries
European_Union=['Austria','Belgium','Bulgaria','Croatia','Cyprus','Czech Republic','Denmark','Estonia',\
                'Finland','France','Germany','Greece','Hungary','Ireland','Italy','Latvia','Lithuania',\
                'Luxembourg','Malta','Netherlands','Poland','Portugal','Romania','Slovakia','Slovenia',\
                'Spain','Sweden','United Kingdom']
#
#--dealing with Serbia and Montenegro
#--remove both countries and only keep Serbia & Montenegro as a whole
emi['Serbia and Montenegro']['co2lu'][2010]=np.sum(emi[country]['co2lu'][2010] for country in ['Serbia','Montenegro'])
del emi['Serbia']
del emi['Montenegro']
World.remove('Serbia')
World.remove('Montenegro')
#--dealing with European Union
for emtype in emtypelist:
  for year in range(1990,2015,5):
    emi['European Union'][emtype][year]=np.sum(emi[country][emtype][year] for country in European_Union)
#
#--now remove individual countries of EU from WorldUnique list
WorldUnique=list(World)   #--has to be listed to be duplicated
for country in European_Union: WorldUnique.remove(country)
#
#--Create country groupings
#--LEA stands for Large Emitters with Absolute target
#--LENA stands for Large Emitters with Non Absolute targets
LEA=['Australia','Brazil','Canada','Japan','Kazakhstan','Russian Federation','Ukraine']
LENA=['Egypt','Indonesia','Iran, Islamic Republic of','Korea, Republic of', 'Malaysia','Mexico','Saudi Arabia','South Africa',\
      'China_Taiwan','Thailand','Turkey','United Arab Emirates']
Transport=['Int. Aviation','Int. Shipping']
#
MainCountries=['United States','European Union','China','India']+LEA+LENA
#
#--Defining remaining countries
WorldOther=list(WorldUnique)
for country in MainCountries+Transport: WorldOther.remove(country)
#
#--Dealing with Other countries
#--Here define with other groupings and rest of rest 
OtherAnnex1=['Andorra','Belarus','Iceland','Liechtenstein','Monaco','New Zealand','Norway','Switzerland']
OtherEmerging=['Chile','Philippines','Viet Nam','Singapore']
#--NDCs of countries from Other Annex 1 and Other Emerging are used directly. 
#--Other Oil countries are treated separately, because their NDCs only inlude some sectors.
OtherOil=['Brunei Darussalam','Bahrain','Kuwait','Oman']
RestOfWorld=list(WorldOther)
for country in OtherAnnex1+OtherEmerging+OtherOil: RestOfWorld.remove(country)
#--countries from Rest of World which NDC can be used directly
RestOfWorldNDC=['Afghanistan','Albania','Angola','Argentina','Azerbaijan','Bangladesh','Barbados','Benin','Bhutan',\
                 'Bosnia and Herzegovina','Botswana','Burkina Faso','Burundi','Cambodia','Cameroon',\
                 'Central African Republic','Chad','Colombia','Comoros','Congo','Congo_the Democratic Republic of the',\
                 'Costa Rica','Djibouti','Dominica','Dominican Republic','Equatorial Guinea','Eritrea','Ethiopia','Gambia',\
                 'Georgia','Ghana','Grenada','Guatemala','Guinea','Haiti','Jamaica','Kenya','Kiribati',\
                 "Korea, Democratic People's Republic of",'Lebanon','Liberia','Macedonia, the former Yugoslav Republic of',\
                 'Madagascar','Maldives','Marshall Islands','Mauritania','Mauritius','Micronesia, Federated States of',\
                 'Moldova, Republic of','Mongolia','Morocco','Namibia','Niger','Nigeria','Pakistan','Paraguay','Peru',\
                 'Saint Kitts and Nevis','Saint Lucia','Saint Vincent and the Grenadines','San Marino',\
                 'Sao Tome and Principe','Senegal','Serbia and Montenegro','Seychelles','Solomon Islands','Tajikistan',\
                 'Tanzania_United Republic of','Togo','Trinidad and Tobago','Tunisia','Uganda','Venezuela','Zambia']
#--other countries from Rest of World
RestOfWorldNoNDC=list(RestOfWorld)
for country in RestOfWorldNDC: RestOfWorldNoNDC.remove(country)
#
WorldOtherGroups=['OtherAnnex1','OtherEmerging','OtherOil','RestOfWorld']
#
RepartCountries=['United States','European Union','LEA','LENA','China','India','WorldOther','Transport']
#
#---The list of countries is included in WorldUnique with 169 countries
#---it can be decomposed as MainCountries+OtherAnnex1+OtherEmerging+OtherOil+RestOfWorld+Transport
#---or MainCountries+OtherAnnex1+OtherEmerging+OtherOil+RestOfWorldNDC+RestOfWorldNoNDC+Transport
#-----------------------------------------------------------------------------------------------------------------
#
#--adding UNFCCC data to emi dictionary
for country in MainCountries: emi[country]['co2luUNFCCC']=emiU[country]
#
#--adding dictionary for GHG emissions excluding LULUCF and total GHG emissions within emi
for country in WorldUnique:
  emiGHGcountry={}
  emitotcountry={}
  for year in range(1990,2015,5):
    emiGHGcountry[year]=0.0
    emitotcountry[year]=0.0
  emi[country]['GHGexclLU']=emiGHGcountry
  emi[country]['GHGtot']=emitotcountry
#
#--Calculating GHG emissions without LULUCF based on EDGAR data
for country in WorldUnique:
  for year in range(1990,2015,5):
    emi[country]['GHGexclLU'][year]=np.sum(emi[country][emtype][year]*gwp[emtype] \
                                           for emtype in emtypelist if emtype<>'co2lu')
#
#--Adding LULUCF emissions, except for the US, Canada and Russia (see below)
#--defining a fraction of the carbon sinks considered as anthropogenic
#--Used for the US, Canada and Russia, already for historical data; 
#--can be modified according to assumption on anthropogenic carbon sinks: here constant median value
alpha_sink_ant=0.5
#
for country in WorldUnique:
  if country not in ['United States','Canada','Russian Federation']:
    for year in range(1990,2015,5):
      emi[country]['GHGtot'][year]=emi[country]['GHGexclLU'][year]+emi[country]['co2lu'][year]
#
  else:
#--Dealing with net-net emissions for the US, Canada and Russia
    for year in range(1990,2015,5):
      emi[country]['GHGtot'][year]=emi[country]['GHGexclLU'][year]+emi[country]['co2lu'][year]+\
                                   alpha_sink_ant*(emi[country]['co2luUNFCCC'][year]-emi[country]['co2lu'][year])
#
#--Adding new categories for groups of countries
typeOther=['GHGexclLU','co2lu','GHGtot']
#
#--Sum up all emissions for past (1990-2010)
for group in WorldOtherGroups+['LEA','LENA','Transport','RestOfWorldNDC','RestOfWorldNoNDC','WorldOther','WorldUnique']:
  emigroup={}
  for type in typeOther:
    emigrouptype={}
    for year in range(1990,2015,5):
      if group<>'WorldUnique':
        emigrouptype[year]=np.sum(emi[country][type][year] for country in eval(group))
      else:
        emigrouptype[year]=np.sum(emi[country][type][year] for country in eval(group))/1.e6
    emigroup[type]=emigrouptype
  emi[group]=emigroup
#
#-------------------------------------------------------------------------------------------------------------------
#
#--READING POPULATION DATA: source IIASA database, SSP scenarios (population identical per SSP for every data source).
#--opening population file
#
fidpop=open(path+'pop_data_SSP.csv','r')
nblineheaderP=2
yearspop=[2005,2010,2020,2030,2040]
nbyearpop=len(yearspop)
#
#--All population values are in millions of people.
#
#--creating empty dictionary
pop={}
#
#--reading through the UN data
nbline=0
for fline in fidpop:
  nbline=nbline+1
  if nbline > nblineheaderP:
    #--reading and splitting line
    data=fline.split(';')
    country=data[0]
    popcountry={}
    #--pop dictionary
    for off,scen in zip(range(len(SSPscenarios)),SSPscenarios):
      popcountryscen={}
      for ind in range(0,nbyearpop):
        year=yearspop[ind]
        if data[1+ind]<>'':
          popcountryscen[year] =float(data[1+ind+off*nbyearpop])
        else:
          if verbose: print 'Attention',country, scen, year, ' missing'
          popcountryscen[year] =0.0
      popcountry[scen]=popcountryscen
    pop[country]=popcountry
#
#--dealing with European Union
for scen in SSPscenarios:
  for year in yearspop:
    pop['European Union'][scen][year]=np.sum(pop[country][scen][year] for country in European_Union)
#
#-------------------------------------------------------------------------------------------------------------
#
#--READING GDP DATA: different growth scenarios related to the 5 Shared Socioeconomic Pathways from the IPCC.
#
GDPsources=['CEPII','OECD','IIASA','PIK']
years=[2005,2010,2020,2030,2040] #--years for which GDP scenario values are provided
yearsann=range(1990, 2017,1)
#--opening GDP file
fidCEPII=open(path+'gdp_data_CEPII_cap.csv','r')
fidOECD=open(path+'gdp_data_OECD.csv','r')
fidIIASA=open(path+'gdp_data_IIASA.csv','r')
fidPIK=open(path+'gdp_data_PIK.csv','r')
fidHist=open(path+'gdp_data_WB_History.csv','r')
fidAnn=open(path+'gdp_data_WB_annual.csv','r')
fid=[fidCEPII,fidOECD,fidIIASA,fidHist]
#
#--reading through World Bank annual historical data from 1990 to 2016
#--BE CAREFUL: values in PPP constant 2011 international $ - need to convert to constant billion USD 2005
WorldGDPann=[]
GDPann={}
nbline=0
#--Inflation values for the USA between 2006 and 2011 are used to convert constant 2011 international $ to constant USD 2005
#--Source World Bank: https://data.worldbank.org/indicator/NY.GDP.DEFL.KD.ZG?locations=US&name_desc=false
inflation={2006: 3.072 , 2007: 2.661 , 2008: 1.962, 2009: 0.759, 2010: 1.221, 2011: 2.065}
for fline in fidAnn:
  nbline+=1
  if nbline>1:
    data=fline.split(';')
    country=data[0]
    WorldGDPann.append(country)
    GDPcountry={}
    #--GDP dictionary
    for ind in range(0,len(yearsann)):
      year=yearsann[ind]
      if data[1+ind]<>''and data[1+ind]<>'\r\n':
        value=float(data[1+ind])
        infl_factor=1.0
        for yr in range(2006,2011+1):
          infl_factor*=(1.0 + inflation[yr]/100.0)
        value/=infl_factor                            #--convert to USD 2005
        value/=1000000000.                              #--convert in billion
      else:
        value=0.0
        if verbose: print 'Attention ',country,' GDP',' value', year, ' missing'
      GDPcountry[year]=value
    GDPann[country]=GDPcountry
#
#--GDP data for decadal years used to calculate the intensity baseline (2005) and the intensity target (2030).
#--Source: OECD SSP scenarios, constant billion USD 2005
#--CEPII SSP scenarios, GDP per capita in 2005 PPP
#--IIASA SSP scenarios + World Bank for 2005 ;
#--PIK SSP scenarios + World Bank for ratios country/region, constant billion USD 2005
#
nblineheaderG=[2,2,2,1]
nbyear=len(years)
#
#--reading through the PIK GDP data, available for 32 regions only.
GDPPIK={}
Regions=[]
nbline=0
for fline in fidPIK:
  nbline=nbline+1
  if nbline > 2:
    #--reading and splitting line
    data=fline.split(';')
    region=data[0]
    GDPregion={}
    #--GDP dictionary
    for off,scen in zip(range(len(SSPscenarios)),SSPscenarios):
      GDPregionscen={}
      for ind in range(0,nbyear):
        year=years[ind]
        if data[1+off*(nbyear)+ind]<>'':
          value=float(data[1+off*(nbyear)+ind])
        else:
          value=0.0
          if verbose: print 'Attention ',region,' PIK GDP',' value', year, ' missing'
        GDPregionscen[year]=value
      GDPregion[scen]=GDPregionscen
    GDPPIK[region]=GDPregion
    Regions.append(region)
#
#--reading through the GDP data
GDPprep={} #--preliminary dictionary: following the structure of the input data
GDP={}     #--final dictionary: same elements, but following the structure of the emissions dictionary
GDPrev={}  #--final dictionary: revision of GDP according to historical GDP
RegionsList=[] #--list of PIK regions corresponding to each country
WorldGDP=[]
yearsHist=[1990,1995,2000,2005]
#
for indS,source in zip(range(0,len(GDPsources)),GDPsources):
  GDPsource={}
  nbline=0
  for fline in fid[indS]:
    nbline=nbline+1
    if nbline > nblineheaderG[indS]:
      #--reading and splitting line
      data=fline.split(';')
      country=data[0]
      GDPcountry={}
      #--GDP dictionary
      if indS<>3:
        for off,scen in zip(range(len(SSPscenarios)),SSPscenarios):
          GDPcountryscen={}
          for ind in range(0,nbyear):
            year=years[ind]
            if data[1+off*nbyear+ind]<>''and data[1+off*nbyear+ind]<>'\r\n':
              value=float(data[1+off*nbyear+ind])
            else:
              value=0.0
              if verbose: print 'Attention ',country,' GDP',' value', year, ' missing'
            GDPcountryscen[year]=value
          GDPcountry[scen]=GDPcountryscen
        GDPsource[country]=GDPcountry
      else: 
        WorldGDP.append(country)
        region=data[1]
        RegionsList.append(region)
        for off,scen in zip(range(len(SSPscenarios)),SSPscenarios):
          GDPcountryscen={}
          for ind in range(0,4):
            year=yearsHist[ind]
            if data[2+ind]<>''and data[2+ind]<>'\r\n':
              valueR=float(data[2+ind])
            else:
              valueR=0.0
              if verbose: print 'Attention ',country, ' World Bank GDP',' value', year, ' missing'
            GDPcountryscen[year]=valueR
          for ind in range(1,nbyear):
            year=years[ind]
            if region<>'':
              #--temporary without the country/region ratio
              value=GDPPIK[region][scen][year]
            else:
              value=0.0
              if verbose: print 'Attention ',country, source, ' GDP',' value', year, ' missing'
            GDPcountryscen[year]=value
          GDPcountry[scen]=GDPcountryscen
          del GDPcountryscen
        GDPsource[country]=GDPcountry
        del GDPcountry
  GDPprep[source]=GDPsource
  del GDPsource
#
#--calculating GDP values for PIK scenarios on country level: partial convergence of emissions per capita per region
convyear=2040. #--convergence year
for scen in SSPscenarios:
  for country,indR in zip(WorldGDP,range(0,len(WorldGDP))):
    nbCR=0.
    for count,indRR in zip(WorldGDP,range(0,len(WorldGDP))):
      if RegionsList[indRR]==RegionsList[indR] and count<>country:
        nbCR+=1
    if nbCR>0:
      if RegionsList[indR]<>'' and  pop[country][scen][2005]<>0. and pop[country][scen][2005]<>'':
        popR=0.
        for count,indRR in zip(WorldGDP,range(0,len(WorldGDP))):
          if RegionsList[indRR]==RegionsList[indR]:
            popR+=pop[count][scen][convyear]
        for ind in range(1,nbyear):
          year=years[ind]
          GDPprep['PIK'][country][scen][year]=GDPprep['PIK'][country][scen][2005]**((convyear-year)/(convyear-2005.))\
                                             *GDPPIK[RegionsList[indR]][scen][convyear]**((year-2005.)/(convyear-2005.))\
                                             *pop[country][scen][year]\
                                             /(pop[country][scen][2005])**((convyear-year)/(convyear-2005.))\
                                             /(popR)**((year-2005.)/(convyear-2005.))
      else:
        for ind in range(1,nbyear):
          year=years[ind]
          GDPprep['PIK'][country][scen][year]=0.
#
for scen in SSPscenarios:
  for country,indR in zip(WorldGDP,range(0,len(WorldGDP))):
    nbCR=0.
    for count,indRR in zip(WorldGDP,range(0,len(WorldGDP))):
      if RegionsList[indRR]==RegionsList[indR] and count<>country:
        nbCR+=1
    if nbCR>0:
      if RegionsList[indR]<>'' and  pop[country][scen][2005]<>0. and pop[country][scen][2005]<>'':
        sumR={}
        for ind in range(0,nbyear):
          year=years[ind]
          sumR[year]=0.
          for countr,indRRR in zip(WorldGDP,range(0,len(WorldGDP))):
            if RegionsList[indR]==RegionsList[indRRR]:
              sumR[year]=sumR[year]+GDPprep['PIK'][countr][scen][year]
          if ind>0:
            GDPprep['PIK'][country][scen][year]=GDPprep['PIK'][country][scen][year]+(GDPprep['PIK'][country][scen][year]\
                                                 -GDPprep['PIK'][country][scen][years[ind-1]])\
                                                /(sumR[year]-sumR[years[ind-1]])\
                                                *(GDPPIK[RegionsList[indR]][scen][year]-sumR[year])
#
#--Inverting the GDP data structure
for source in GDPsources:
  GDPsource={}
  for scen in SSPscenarios:
    GDPscen={}
    for country in MainCountries+European_Union+WorldOther:
      GDPscencountry={}
      for year in years:
        GDPscencountry[year]=GDPprep[source][country][scen][year]
      GDPscen[country]=GDPscencountry
      del GDPscencountry
    GDPsource[scen]=GDPscen
    del GDPscen
  GDP[source]=GDPsource
  del GDPsource
#
#--Converting GDP per cap from CEPII into emissions with population data
for scen in SSPscenarios:
  for year in years:
    for country in MainCountries+European_Union+WorldOther:
      GDP['CEPII'][scen][country][year]*=pop[country][scen][year]/1.e3
#
#--Summing GDP data for countries in European Union
for source in GDPsources:
  for scen in SSPscenarios:
    for year in years:
      GDP[source][scen]['European Union'][year]=np.sum(GDP[source][scen][country][year] for country in European_Union)
#
#--Building a time series of GDP for countries with intensity targets
#--which is constrained by past values of GDP data if revise_GDP is True
#--and only fitted to the scenarios if revise_GDP is False
GDPyears=[2005,2010,2020,2030,2040]
for source in GDPsources:
  GDPsource={}
  for scen in SSPscenarios:
    GDPscen={}
    for country in WorldGDPann:
      GDPcountry={}
      GROWTHcountry={}
      #--We first check if the actual (historical) GDP exists for that country
      GDPann_exists=True
      for year in range(2005,2016):
        if GDPann[country][year] < 1.e-10: GDPann_exists=False
      if not GDPann_exists  and verbose: print 'Be careful, no actual GDP data exist for ', country
      #--We then check if the SSP GDP exists for that country
      GDPscen_exists=True
      for year in GDPyears:
        if GDP[source][scen][country][year] < 1.e-10: GDPscen_exists=False
      if not GDPscen_exists and verbose: print 'Be careful, no scenario GDP data exist for ', country
      #--if the SSP projection data do not exist for a country we just assign zero to GDP values
      if not GDPscen_exists:
        for year in range(2005,2041):
          GDPcountry[year]=0.0
      #--the SSP projection data do exist
      else:
        #--the actual GDP exists so we constrain the scenario data
        if revise_GDP and GDPann_exists: 
          for year in range(2005,2016):
            #--We use the actual GDP and compute the growth
            GDPcountry[year]=GDPann[country][year]
            #--Compute growth
            if year > 2005: 
              GROWTHcountry[year]=GDPcountry[year]/GDPcountry[year-1]-1.0
          #--fit GDP data from scenarios
          GDPdata=[]
          for yr in GDPyears:
            GDPdata.append(GDP[source][scen][country][yr])
          GDPfit=interpolate.interp1d(GDPyears, GDPdata, kind='cubic')
          #--now building the revised GDP data
          for year in range(2016,2041):
            #--growth from the scenario
            growth_scen = GDPfit(year)/GDPfit(year-1) - 1.0
            #--weights
            if year <= 2020: 
               x1, x2, x3 = 1.0, 1.0, 1.0
            elif year <= 2030:
               x1 = 1.0 + 2.0*float(year-2020)/float(2030-2020)
               x2 = 1.0 - 1.0*float(year-2020)/float(2030-2020)
               x3 = 1.0 - 1.0*float(year-2020)/float(2030-2020)
            else:
               x1, x2, x3 = 3.0, 0.0, 0.0 
            #--computing weighted growth
            GROWTHcountry[year] = (x1*growth_scen+x2*GROWTHcountry[year-1]+x3*GROWTHcountry[year-2])/(x1+x2+x3)
            #--and the GDP
            GDPcountry[year]=GDPcountry[year-1]*(1.+GROWTHcountry[year])
        #
        #--the actual GDP does not exists or there is no revision so we simply interpolate the scenario data
        else:
          GDPdata=[]
          for yr in GDPyears:
            GDPdata.append(GDP[source][scen][country][yr])
          GDPfit=interpolate.interp1d(GDPyears, GDPdata, kind='cubic')
          for year in range(2005,2041):
            GDPcountry[year]=GDPfit(year)
  #--packing up
      GDPscen[country]=GDPcountry
      del GDPcountry
    GDPsource[scen]=GDPscen
    del GDPscen
  GDPrev[source]=GDPsource
  del GDPsource
#
#-----------------------------------------------------------------------------------------------------------------------
#
#--Gathering informations from the NDCs
#
#--NDCs expressed for target years other than 2030 are directly converted for 2030.
#--Summarizing most ambitious versions of NDCs in terms of reduction percentage or absolute values
NDC_high={'United States':-0.28,'European Union':-0.40,'China':-0.65,'Canada':-0.30,'Mexico':-0.36,\
           'Brazil':-0.43,'South Africa':398000,'Russian Federation':-0.30,'Japan':-0.25,\
           'Korea, Republic of':-0.37,'Australia':-0.28,'India':-0.35,'Indonesia':-0.41,\
           'United Arab Emirates':-0.40,'Egypt':-0.40,'Iran, Islamic Republic of':-0.40,'Saudi Arabia':-0.40,\
           'Kazakhstan':-0.25,'Malaysia':-0.45,'China_Taiwan':-0.50,'Turkey':-0.21,'Thailand':-0.25,\
           'Ukraine':-0.40,'Int. Aviation':906000,'Int. Shipping':940000,'Andorra':-0.37,'Belarus':-0.28,\
           'Iceland':-0.40,'Liechtenstein':-0.40,'Monaco':-0.50,'New Zealand':-0.30,'Norway':-0.40,'Switzerland':-0.50,\
           'Chile':-0.45,'Philippines':-0.45,'Singapore':-0.36,'Viet Nam':-0.25,'OtherOil':-0.40,\
           'Afghanistan':-0.136,'Albania':-0.115,'Angola':-0.5,'Argentina':-0.37,'Azerbaijan':-0.35,\
           'Bangladesh':-0.15,'Barbados':-0.44,'Benin':-0.214,'Bhutan':6300,'Bosnia and Herzegovina':-0.23,\
           'Botswana':-0.15,'Burkina Faso':-0.3695,'Burundi':-0.20,'Cambodia':-0.27,'Cameroon':-0.32,\
           'Central African Republic':-0.05,'Chad':-0.71,'Colombia':-0.30,'Comoros':-0.84,'Congo':-0.51,\
           'Congo_the Democratic Republic of the':-0.17,'Costa Rica':9370,\
           'Djibouti':-0.60,'Dominica':-0.447,'Dominican Republic':-0.25,'Equatorial Guinea':-0.2,\
           'Eritrea':-0.39,'Ethiopia':-0.64,'Gambia':-0.82,'Georgia':-0.25,'Ghana':-0.45,'Grenada':-0.40,\
           'Guatemala':-0.112,'Guinea':-0.13,'Haiti':-0.26,'Jamaica':-0.10,\
           'Kenya':-0.3,'Kiribati':-0.618,"Korea, Democratic People's Republic of":-0.4025,\
           'Lebanon':-0.3,'Liberia':-0.15,'Macedonia, the former Yugoslav Republic of':-0.36,'Madagascar':-0.14,\
           'Maldives':-0.24,'Marshall Islands':-0.45,'Mauritania':-0.1962,'Mauritius':-0.30,\
           'Micronesia, Federated States of':-0.35,'Moldova, Republic of':-0.78,'Mongolia':-0.14,'Morocco':-0.42,'Namibia':-0.89,\
           'Niger':-0.346,'Nigeria':-0.45,'Pakistan':-0.2,'Paraguay':-0.2,'Peru':-0.3,'Saint Kitts and Nevis':-0.35,\
           'Saint Lucia':-0.23,'Saint Vincent and the Grenadines':-0.22, 'San Marino':-0.2,'Sao Tome and Principe':-0.24,\
           'Senegal':-0.21,'Serbia and Montenegro':-0.098,'Seychelles':-0.29,'Solomon Islands':-0.45,'Tajikistan':-0.35,\
           'Tanzania_United Republic of':-0.20,'Togo':-0.3114,'Trinidad and Tobago':-0.15,'Tunisia':-0.41,\
           'Uganda':-0.22,'Venezuela':-0.2,'Zambia':-0.88,'RestOfWorldNoNDC':-0.45}
#
#--Summarizing least ambitious versions of NDCs in terms of reduction percentage or absolute values
NDC_low ={'United States':-0.26,'European Union':-0.40,'China':-0.60,'Canada':-0.30,'Mexico':-0.22,\
           'Brazil':-0.43,'South Africa':614000,'Russian Federation':-0.25,'Japan':-0.25,\
           'Korea, Republic of':-0.37,'Australia':-0.26,'India':-0.33,'Indonesia':-0.29,\
           'United Arab Emirates':-0.30,'Egypt':-0.30,'Iran, Islamic Republic of':-0.30,'Saudi Arabia':-0.30,\
           'Kazakhstan':-0.15,'Malaysia':-0.35,'China_Taiwan':-0.50,'Turkey':-0.21,'Thailand':-0.20,\
           'Ukraine':-0.40,'Int. Aviation':1200000,'Int. Shipping':1200000,'Andorra':-0.37,'Belarus':-0.28,\
           'Iceland':-0.40,'Liechtenstein':-0.40,'Monaco':-0.50,'New Zealand':-0.30,'Norway':-0.40,'Switzerland':-0.50,
           'Chile':-0.30,'Philippines':-0.30,'Singapore':-0.36,'Viet Nam':-0.08,'OtherOil':-0.30,\
           'Afghanistan':0.,'Albania':-0.115,'Angola':-0.35,'Argentina':-0.18,'Azerbaijan':-0.35,\
           'Bangladesh':-0.05,'Barbados':-0.44,'Benin':-0.035,'Bhutan':6300,'Bosnia and Herzegovina':-0.02,\
           'Botswana':-0.15,'Burkina Faso':-0.066,'Burundi':-0.03,'Cambodia':-0.27,'Cameroon':0.,\
           'Central African Republic':0.,'Chad':-0.182,'Colombia':-0.20,'Comoros':-0.84,'Congo':0.,\
           'Congo_the Democratic Republic of the':-0.,'Costa Rica':9370,\
           'Djibouti':-0.40,'Dominica':-0.447,'Dominican Republic':-0.,'Equatorial Guinea':-0.2,\
           'Eritrea':-0.,'Ethiopia':-0.,'Gambia':-0.13,'Georgia':-0.15,'Ghana':-0.15,'Grenada':-0.40,\
           'Guatemala':-0.112,'Guinea':-0.,'Haiti':-0.05,'Jamaica':-0.08,\
           'Kenya':-0.,'Kiribati':-0.128,"Korea, Democratic People's Republic of":-0.08,\
           'Lebanon':-0.15,'Liberia':-0.15,'Macedonia, the former Yugoslav Republic of':-0.3,'Madagascar':-0.14,\
           'Maldives':-0.10,'Marshall Islands':-0.45,'Mauritania':-0.0268,'Mauritius':-0.,\
           'Micronesia, Federated States of':-0.28,'Moldova, Republic of':-0.65,'Mongolia':-0.14,'Morocco':-0.17,'Namibia':-0.,\
           'Niger':-0.035,'Nigeria':-0.20,'Pakistan':-0.,'Paraguay':-0.1,'Peru':-0.2,'Saint Kitts and Nevis':-0.35,\
           'Saint Lucia':-0.23,'Saint Vincent and the Grenadines':-0.22, 'San Marino':-0.2,'Sao Tome and Principe':-0.,\
           'Senegal':-0.05,'Serbia and Montenegro':-0.098,'Seychelles':-0.,'Solomon Islands':-0.30,'Tajikistan':-0.10,\
           'Tanzania_United Republic of':-0.10,'Togo':-0.1114,'Trinidad and Tobago':-0.0025,'Tunisia':-0.13,\
           'Uganda':-0.22,'Venezuela':-0.2,'Zambia':-0.35,'RestOfWorldNoNDC':-0.30}
#
#--Explaining which type of target is used in the NDC: absolute, in carbon intensity of GDP \
#--compared to a BAU or with an emission value directly given in the NDC.
NDC_type={'United States':'absolute','European Union':'absolute','China':'intensity','Canada':'absolute','Mexico':'BAU',\
           'Brazil':'absolute','South Africa':'value','Russian Federation':'absolute','Japan':'absolute',\
           'Korea, Republic of':'BAU','Australia':'absolute','India':'intensity','Indonesia':'BAU',\
           'United Arab Emirates':'intensity','Egypt':'intensity','Iran, Islamic Republic of':'intensity',\
           'Saudi Arabia':'intensity','Kazakhstan':'absolute','Malaysia':'intensity','China_Taiwan':'BAU',\
           'Turkey':'BAU','Thailand':'BAU','Ukraine':'absolute','Int. Aviation':'value','Int. Shipping':'value',\
           'Andorra':'BAU','Belarus':'absolute','Iceland':'absolute','Liechtenstein':'absolute','Monaco':'absolute',\
           'New Zealand':'absolute','Norway':'absolute','Switzerland':'absolute',
           'Chile':'intensity','Philippines':'intensity','Singapore':'intensity','Viet Nam':'BAU','OtherOil':'intensity',\
           'Afghanistan':'BAU','Albania':'BAU','Angola':'BAU','Argentina':'BAU','Azerbaijan':'absolute',\
           'Bangladesh':'BAU','Barbados':'BAU','Benin':'BAU','Bhutan':'value','Bosnia and Herzegovina':'BAU',\
           'Botswana':'absolute','Burkina Faso':'BAU','Burundi':'BAU','Cambodia':'BAU','Cameroon':'BAU',\
           'Central African Republic':'BAU','Chad':'BAU','Colombia':'BAU','Comoros':'BAU','Congo':'BAU',\
           'Congo_the Democratic Republic of the':'BAU','Costa Rica':'value',\
           'Djibouti':'BAU','Dominica':'absolute','Dominican Republic':'absolute',\
           'Equatorial Guinea':'absolute','Eritrea':'BAU','Ethiopia':'BAU','Gambia':'absolute','Georgia':'BAU',\
           'Ghana':'BAU','Grenada':'absolute','Guatemala':'BAU','Guinea':'absolute','Haiti':'BAU','Jamaica':'BAU',\
           'Kenya':'BAU','Kiribati':'BAU',"Korea, Democratic People's Republic of":'BAU',\
           'Lebanon':'BAU','Liberia':'BAU','Macedonia, the former Yugoslav Republic of':'BAU','Madagascar':'BAU',\
           'Maldives':'BAU','Marshall Islands':'absolute','Mauritania':'BAU','Mauritius':'BAU',\
           'Micronesia, Federated States of':'absolute','Moldova, Republic of':'absolute','Mongolia':'BAU','Morocco':'BAU','Namibia':'BAU',\
           'Niger':'BAU','Nigeria':'BAU','Pakistan':'BAU','Paraguay':'BAU','Peru':'BAU','Saint Kitts and Nevis':'BAU',\
           'Saint Lucia':'BAU','Saint Vincent and the Grenadines':'BAU','San Marino':'absolute','Sao Tome and Principe':'BAU',\
           'Senegal':'BAU','Serbia and Montenegro':'absolute','Seychelles':'BAU','Solomon Islands':'BAU',\
           'Tajikistan':'absolute','Tanzania_United Republic of':'BAU','Togo':'BAU','Trinidad and Tobago':'BAU',\
           'Tunisia':'intensity','Uganda':'BAU','Venezuela':'BAU','Zambia':'BAU','RestOfWorldNoNDC':'intensity'}
#
#--We choose to define the remainings of those groups of countries' targets in terms of intensity, as for China and India.
#
#--Summarizing base year or BAU value in 2030.
#Countries giving their NDC in emission values have their baseline set to 0.
Baseline ={'United States':2005,'European Union':1990,'China':2005,'Canada':2005,'Mexico':973000,\
           'Brazil':2005,'South Africa':0,'Russian Federation':1990,'Japan':2005,\
           'Korea, Republic of':850600,'Australia':2005,'India':2005,'Indonesia':2869000,\
           'United Arab Emirates':2005,'Egypt':2005,'Iran, Islamic Republic of':2005,\
           'Saudi Arabia':2005,'Kazakhstan':1990,'Malaysia':2005,'China_Taiwan':428000,\
           'Turkey':1175000,'Thailand':555000,'Ukraine':1990,'Int. Aviation':0,
           'Int. Shipping':0,'Andorra':530.55,'Belarus':1990,'Iceland':1990,'Liechtenstein':1990,\
           'Monaco':1990,'New Zealand':2005,'Norway':1990,'Switzerland':1990,\
           'Chile':2005,'Philippines':2005,'Singapore':2005,'Viet Nam':787400,'OtherOil':2005,\
           'Afghanistan':48900,'Albania':6157,'Angola':193250,'Argentina':592000,'Azerbaijan':1990,\
           'Bangladesh':234000,'Barbados':2500,'Benin':22500,'Bhutan':0,'Bosnia and Herzegovina':43000,\
           'Botswana':2010,'Burkina Faso':118300,'Burundi':75000,'Cambodia':11600,'Cameroon':104000,\
           'Central African Republic':100000,'Chad':28700,'Colombia':330000,'Comoros':523,'Congo':27450,\
           'Congo_the Democratic Republic of the':430000,'Costa Rica':0,\
           'Djibouti':4500,'Dominica':2010,'Dominican Republic':2010,\
           'Equatorial Guinea':2010,'Eritrea':7900,'Ethiopia':398000,'Gambia':2010,'Georgia':38420,\
           'Ghana':73950,'Grenada':2010,'Guatemala':54000,'Guinea':1990,'Haiti':45000,'Jamaica':13750,\
           'Kenya':143000,'Kiribati':300,"Korea, Democratic People's Republic of":187730,\
           'Lebanon':43000,'Liberia':5000,'Macedonia, the former Yugoslav Republic of':18000,'Madagascar':214300,\
           'Maldives':3000,'Marshall Islands':2010,'Mauritania':18840,'Mauritius':7000,\
           'Micronesia, Federated States of':2000,'Moldova, Republic of':1990,'Mongolia':52000,'Morocco':170800,'Namibia':22500,\
           'Niger':96468,'Nigeria':900000,'Pakistan':1603000,'Paraguay':416000,'Peru':298000,'Saint Kitts and Nevis':836,\
           'Saint Lucia':817,'Saint Vincent and the Grenadines':600,'San Marino':2005,'Sao Tome and Principe':237,\
           'Senegal':36000,'Serbia and Montenegro':1990,'Seychelles':650,'Solomon Islands':70,\
           'Tajikistan':1990,'Tanzania_United Republic of':153000,'Togo':39000,'Trinidad and Tobago':680000,\
           'Tunisia':2010,'Uganda':77000,'Venezuela':340000,'Zambia':800000,'RestOfWorldNoNDC':2005}
#
#--Mentioning if the target comes directly from the NDC, if supplemental hypothsesis were made
#  or if just a best guess was used.
#
Source={'United States':'suppl hypo','European Union':'direct NDC','China':'suppl hypo',\
        'Canada':'direct NDC','Mexico':'direct NDC','Brazil':'direct NDC','South Africa':'suppl hypo',\
        'Russian Federation':'direct NDC','Japan':'direct NDC','Korea, Republic of':'direct NDC',\
        'Australia':'direct NDC','India':'suppl hypo','Indonesia':'direct NDC',\
        'United Arab Emirates':'best guess','Egypt':'best guess','Iran, Islamic Republic of':'best guess',\
        'Saudi Arabia':'best guess','Kazakhstan':'direct NDC','Malaysia':'suppl hypo',
        'China_Taiwan':'direct NDC','Turkey':'direct NDC','Thailand':'direct NDC','Ukraine':'direct NDC',\
        'Int. Aviation':'best guess','Int. Shipping':'best guess',\
        'Andorra':'direct NDC','Belarus':'direct NDC','Iceland':'direct NDC','Liechtenstein':'direct NDC',\
        'Monaco':'direct NDC','New Zealand':'direct NDC','Norway':'direct NDC','Switzerland':'direct NDC',\
        'Chile':'suppl hypo','Philippines':'best guess','Singapore':'direct NDC','Viet Nam':'suppl hypo',\
        'OtherOil':'best guess',\
        'Afghanistan':'direct NDC','Albania':'direct NDC','Angola':'direct NDC','Argentina':'direct NDC',\
        'Azerbaijan':'direct NDC','Bangladesh':'direct NDC','Barbados':'direct NDC','Benin':'direct NDC',\
        'Bhutan':'direct NDC','Bosnia and Herzegovina':'direct NDC','Botswana':'direct NDC','Burkina Faso':'direct NDC',\
        'Burundi':'direct NDC','Cambodia':'direct NDC','Cameroon':'suppl hypo','Central African Republic':'direct NDC',\
        'Chad':'direct NDC','Colombia':'direct NDC','Comoros':'direct NDC','Congo':'direct NDC',\
        'Congo_the Democratic Republic of the':'direct NDC','Costa Rica':'direct NDC','Djibouti':'direct NDC',\
        'Dominica':'suppl hypo','Dominican Republic':'direct NDC','Equatorial Guinea':'direct NDC','Eritrea':'direct NDC',\
        'Ethiopia':'direct NDC','Gambia':'direct NDC','Georgia':'direct NDC','Ghana':'direct NDC','Grenada':'direct NDC',\
        'Guatemala':'direct NDC','Guinea':'suppl hypo','Haiti':'direct NDC','Jamaica':'direct NDC','Kenya':'direct NDC',\
        'Kiribati':'direct NDC',"Korea, Democratic People's Republic of":'direct NDC','Lebanon':'direct NDC',\
        'Liberia':'direct NDC','Macedonia, the former Yugoslav Republic of':'direct NDC','Madagascar':'direct NDC',\
        'Maldives':'direct NDC','Marshall Islands':'direct NDC','Mauritania':'direct NDC','Mauritius':'direct NDC',\
        'Micronesia, Federated States of':'direct NDC','Moldova, Republic of':'suppl hypo','Mongolia':'direct NDC',\
        'Morocco':'direct NDC','Namibia':'direct NDC','Niger':'direct NDC','Nigeria':'direct NDC','Pakistan':'direct NDC',\
        'Paraguay':'direct NDC','Peru':'direct NDC','Saint Kitts and Nevis':'direct NDC',\
        'Saint Lucia':'direct NDC','Saint Vincent and the Grenadines':'suppl hypo','San Marino':'direct NDC',\
        'Sao Tome and Principe':'direct NDC','Senegal':'direct NDC','Serbia and Montenegro':'suppl hypo',\
        'Seychelles':'direct NDC','Solomon Islands':'direct NDC','Tajikistan':'direct NDC',\
        'Tanzania_United Republic of':'direct NDC','Togo':'direct NDC','Trinidad and Tobago':'direct NDC',\
        'Tunisia':'suppl hypo','Uganda':'direct NDC','Venezuela':'direct NDC','Zambia':'direct NDC',\
        'RestOfWorldNoNDC':'best guess'}
# 
#--Target for the Philippines, part of Other Emerging, is based on the Chilean NDC: no BAU value given. \
#--Targets for Other Oil and Rest of World No NDC groups are not based on NDCs, mostly containing only sectoral targets.
#--Targets for international aviation and shipping are based on external sources
#
#--Gathering all infos on NDCs for countries in one dictionary
NDC={}
for country in MainCountries+Transport+OtherAnnex1+OtherEmerging+RestOfWorldNDC+['OtherOil','RestOfWorldNoNDC']:
  NDC[country]={'high':NDC_high[country], 'low':NDC_low[country], \
                 'type':NDC_type[country], 'baseline':Baseline[country], 'source':Source[country]}
#
#-------------------------------------------------------------------------------------------------------------------
#--Dealing with the emission peak constraint for China
#
#--Estimating GDP values per year with polynomial fit
yearsIT=range(2011,2031,1) #--period over which we compute the interpolation of the carbon intensity for China
#
GDPyears=[2005,2010,2020,2030,2040]
GDPChinagrowth={} #--Chinese GDP growth dictionary
for source in GDPsources:
  GDPChinagrowth[source]={}
  for scen in SSPscenarios:
    GDPChinagrowth[source][scen]={}
    for year in yearsIT:
      GDPChinagrowth[source][scen][year]=GDPrev[source][scen]['China'][year]/GDPrev[source][scen]['China'][year-1]-1.0
#
#--Energy Intensity for China
IT={}
for source in GDPsources:
  IT[source]={}
  for scen in SSPscenarios:
    IT[source][scen]={}
    for year in [2005,2010]:
      IT[source][scen][year]=emi['China']['co2'][year]/GDPrev[source][scen]['China'][year]
#
def ITfunction(source,scen,irt,year):
#
  if year in range(2006,2030):
    IT_2030=(1.+irt)*IT[source][scen][2005]
    IT_tcam=IT[source][scen][2005]*((IT_2030/IT[source][scen][2005])**(1./25.))**(year-2005)
    IT_tlin=IT[source][scen][2005]+(year-2005.)/(2030.-2005.)*(IT_2030-IT[source][scen][2005])
  else:
    print 'Error in year'
    quit()
  return IT_tcam, IT_tlin
#
#-------------------------------------------------------------------------------------------------------------------------------------
#
#--Parameters considered in probabilistic terms
#
#--Parameters related to the ambition level of NDCs
#--NDCs high and low assumptions as defined above, both unconditional and conditional statements in NDCs, and our assumptions
#
#--Parameters related to Land Use emissions consideration and evolution
#--The Indian NDC indicates a target of creating an additional carbon sink of 2.5 to 3 billion tCO2eq through
#--additional forest and tree cover by 2030. We therefore divide this target over the next fifteen years, and
#--add it to the carbon sink of 2010. 
#
emi_India_co2lu={}
emi_India_co2lu[2030]={'high':emi['India']['co2lu'][2010]-3000000./15., 'low':emi['India']['co2lu'][2010]-2500000./15.}
#
#--Hypothesis for LULUCF emissions reduction in 2030 compared to 2010 for country group Rest of World
ROWredLUhigh=0.0
ROWredLUlow =0.5
#
#--Parameters related to non-CO2 emissions from China
#--Ratio CO2eq/CO2 for China: extrapolation from 1990-2010, following a fit in log(CO2 emissions)
ratioChina={}
ratioC=[]
co2C=[]
for year in range(1990,2015,5):
  co2exclLUChina=emi['China']['co2'][year]
  co2C.append(co2exclLUChina)
  ratioC.append(emi['China']['GHGexclLU'][year]/co2exclLUChina)
ratioCfit=np.polyfit(np.log(co2C),ratioC,1,full=True)
#--considering a range of [-sigma,+sigma] around the extrapolated value of the ratio,\
#--with sigma given by the Python least-squared fit
emi_China_nonCO2={'high':(ratioCfit[1][0]/len(range(1990,2015,5)))**(0.5),\
                  'low':-(ratioCfit[1][0]/len(range(1990,2015,5)))**(0.5)}
#
#--Parameters related to economic growth scenarios
#--GDP scenarios sources: OECD, CEPII, IIASA, PIK. Parameter already defined above when reading GDP data
#--GDP scenarios SSP types: SSP1, SSP2, SSP3, SSP4, SSP5. Parameter already defined above when reading GDP data. \
#
#-----------------------------------------------------------------------------------------------------------------------------
#

print("Data loaded")
print("Starting computation. Please wait...")

#--Initialize outputs
#
#--Lists of the emissions results of each run of the Monte-Carlo
List_WorldEmi={}   #--dictionary to store total world emissions from Monte-Carlo
List_Emi={}        #--dictionary to store emissions from individual countries or groupings
List_ChinaEmi={}   #--dictionary to store emissions from China
List_ChinaIRT={}   #--dictionary to store Intensity Reduction for China
#
#--Building the dictionary structure
for source in GDPsources:
  List_WorldEmi[source]={}
  List_Emi[source]={}
  List_ChinaEmi[source]={}
  List_ChinaIRT[source]={}
  for scen in SSPscenarios:
    List_WorldEmi[source][scen]={}
    List_Emi[source][scen]={}
    List_ChinaEmi[source][scen]={}
    List_ChinaIRT[source][scen]=[]
    for country in WorldUnique+WorldOtherGroups+['LEA','LENA','Transport','WorldOther']:
      List_Emi[source][scen][country]={}
      for emtype in typeOther:
        List_Emi[source][scen][country][emtype]={}
        for year in [2030]:
          List_Emi[source][scen][country][emtype][year]=[]
    for emtype in ['co2']+typeOther:
      List_ChinaEmi[source][scen][emtype]={}
      for year in [2030]:
        List_ChinaEmi[source][scen][emtype][year]={}
        for c in ['nopeak','peak']:
          List_ChinaEmi[source][scen][emtype][year][c]=[]
    for emtype in typeOther:
      List_WorldEmi[source][scen][emtype]={}
      for year in [2030]:
        List_WorldEmi[source][scen][emtype][year]=[]
#
#--Complete emissions dictionary with 2030 values
for country in WorldUnique+Transport+['RestOfWorldNDC','RestOfWorldNoNDC']+WorldOtherGroups:
  for emtype in typeOther:
    emi[country][emtype][2030]={}
    for source in GDPsources:
      emi[country][emtype][2030][source]={}
      for scen in SSPscenarios:
        if country=='China': #--forcing the peak by 2030 or not
          emi[country][emtype][2030][source][scen]={}
          for c in ['nopeak','peak']:
            emi[country][emtype][2030][source][scen][c]=0.
        else:
          emi[country][emtype][2030][source][scen]=0.
#
emi['China']['co2'][2030]={}  #--adding the details of CO2 emissions without LULUCF for China
for source in GDPsources:
  emi['China']['co2'][2030][source]={}
  for scen in SSPscenarios:
    emi['China']['co2'][2030][source][scen]={}
    for c in ['nopeak','peak']:
      emi['China']['co2'][2030][source][scen][c]=0.
#
#-------------------------------------------------------------------------------------------------
#
#--Run Monte-Carlo with Ncase members
#
Ncase=50000
ROWredLUrandom=[]
emiIndiaLUrandom=[]
emiChinanonco2random=[]
NDCrandom=[]
ITrandom=[]
#
#--main loop over cases of Monte Carlo simulations
#
for case in range(Ncase):
  #
  #--pick up random values for uncertain parameters
  ROWredLUrandom.append(random())
  ROWredLU=ROWredLUlow+ROWredLUrandom[case]*(ROWredLUhigh-ROWredLUlow)
  #
  emiIndiaLUrandom.append(random())
  #
  emiChinanonco2random.append(random())
  emiChinanonco2=emi_China_nonCO2['low']+emiChinanonco2random[case]*(emi_China_nonCO2['high']-emi_China_nonCO2['low'])
  #
  ITrandom.append(random())
  #
  NDCrandom.append(random())
  for country in MainCountries+Transport+OtherAnnex1+OtherEmerging+RestOfWorldNDC+['OtherOil','RestOfWorldNoNDC']:
    NDC[country]['value']=NDC[country]['low']+NDCrandom[case]*(NDC[country]['high']-NDC[country]['low'])
  #
  #------------------------------------------------------------------------------------------------------
  #
  #--Predict 2030 emissions
  #
  for source in GDPsources:
    #
    for scen in SSPscenarios:
      #
      #--Dealing with biggest countries separately
      #
      for country in MainCountries+Transport:
        #
        #--absolute NDC
        #
        if NDC[country]['type']=='absolute':
          #
          #--implementing net-net approach for USA, Canada and Russia and assuming that national anthropic carbon sinks stay 
          #--constant compared to average over 2000-2010.
          #
          if country in ['United States','Russian Federation']: 
            meanFAO   =np.mean([emi[country]['co2lu'][year] for year in [2000, 2005, 2010]])
            meanUNFCCC=np.mean([emi[country]['co2luUNFCCC'][year] for year in [2000,2005,2010]])
            if country=='United States':#--target given for 2025, linear extrapolation to 2030
              emi[country]['GHGtot'][2030][source][scen]=(emi[country]['GHGexclLU'][NDC[country]['baseline']]+\
                                                       emi[country]['co2luUNFCCC'][NDC[country]['baseline']])* \
                                                       (1.+NDC[country]['value']*(2030-NDC[country]['baseline'])/\
                                                       (2025-NDC[country]['baseline']))-(1-alpha_sink_ant)*(meanUNFCCC-meanFAO)
            else:
              emi[country]['GHGtot'][2030][source][scen]=(emi[country]['GHGexclLU'][NDC[country]['baseline']]+\
                                                       emi[country]['co2luUNFCCC'][NDC[country]['baseline']])* \
                                                       (1.+NDC[country]['value'])-(1-alpha_sink_ant)*(meanUNFCCC-meanFAO)
          #
          #--implementing Canadian net-net approach and assuming that national anthropic carbon sinks stay 
          #--constant compared to average over 2000-2010.
          #
          elif country=='Canada': 
            meanFAO   =np.mean([emi[country]['co2lu'][year] for year in [2000, 2005, 2010]])
            meanUNFCCC=np.mean([emi[country]['co2luUNFCCC'][year] for year in [2000,2005,2010]])
            emi[country]['GHGtot'][2030][source][scen]=emi[country]['GHGexclLU'][NDC[country]['baseline']]*\
                                                      (1.+NDC[country]['value'])+\
                                                       emi[country]['co2luUNFCCC'][NDC[country]['baseline']]\
                                                      -(1-alpha_sink_ant)*(meanUNFCCC-meanFAO)
          #
          #--Micronesia: target year is 2025, linear extrapolation to 2030
          elif country=='Micronesia':
            emi[country]['GHGtot'][2030][source][scen]=emi[country]['GHGtot'][NDC[country]['baseline']]*\
                                                      (1.+NDC[country]['value'])*(2030-NDC[country]['baseline'])/\
                                                      (2025-NDC[country]['baseline'])
          #
          #--other absolute NDC
          # 
          else:
            emi[country]['GHGtot'][2030][source][scen]=emi[country]['GHGtot'][NDC[country]['baseline']]*(1.+NDC[country]['value'])
        #
        #--BaU NDC
        #
        elif NDC[country]['type']=='BAU':
          emi[country]['GHGtot'][2030][source][scen]=NDC[country]['baseline']*(1.+NDC[country]['value'])
        #
        #--value NDC
        #
        elif NDC[country]['type']=='value':
          emi[country]['GHGtot'][2030][source][scen]=NDC[country]['value']
        #
        #--intensity NDC
        #
        elif NDC[country]['type']=='intensity':
          #--for China, India and other countries with intensity targets, the intensity target is applied to GHG emissions without LULUCF.
          #--Except when otherwise mentioned: Malaysia.
          #--India and main oil exporting countries:
          #
          #--Malaysia: we assume that emissions from non-forest land are included in the target.
          if country=='Malaysia':
            emi[country]['GHGtot'][2030][source][scen]=emi[country]['GHGtot'][NDC[country]['baseline']]/         \
                                                          GDPrev[source][scen][country][NDC[country]['baseline']]*  \
                                                          GDPrev[source][scen][country][2030]*(1.+NDC[country]['value'])
          elif country<>'China':
            emi[country]['GHGexclLU'][2030][source][scen]=emi[country]['GHGexclLU'][NDC[country]['baseline']]/    \
                                                          GDPrev[source][scen][country][NDC[country]['baseline']]*   \
                                                          GDPrev[source][scen][country][2030]*(1.+NDC[country]['value'])
            if country=='India': #--based on Indian NDC
              emi[country]['co2lu'][2030][source][scen]=emi_India_co2lu[2030]['low']+  \
                                                      emiIndiaLUrandom[case]*(emi_India_co2lu[2030]['high']-emi_India_co2lu[2030]['low'])
            else: #--We assume LULUCF emissions in 2030 at 2010 level
              emi[country]['co2lu'][2030][source][scen]=emi[country]['co2lu'][2010]
            #
            emi[country]['GHGtot'][2030][source][scen]=emi[country]['GHGexclLU'][2030][source][scen]\
                                                      +emi[country]['co2lu'][2030][source][scen]
      #
      #--China is now treated separately
      #--for China, target is applied to CO2 emissions only, according to the NDC.
      #--first, we compute the minimum intensity reduction target that would imply a peak in Chinese CO2 emissions by 2030:
      #
      #--copy irt value from ambition level
      irt=NDC['China']['value']
      #--estimate corresponding IT value for 2030
      IT_2030=(1.+irt)*IT[source][scen][2005]
      #--calculate the implied intensity value for 2029 as the weighted average between linear and exponential interpolations
      IT_tcam, IT_tlin = ITfunction(source,scen,irt,2029)
      IT_2029 = ITrandom[case]*IT_tcam+(1.-ITrandom[case])*IT_tlin
      #
      #--increase intensity reduction until it is larger than GDPgrowth in 2030
      while -(IT_2030-IT_2029)/IT_2029 < GDPChinagrowth[source][scen][2030] :
        irt=irt-0.002
        IT_2030=(1.+irt)*IT[source][scen][2005]
        #--calculate again the implied intensity value for 2029 as the weighted average between linear and exponential interpolations
        IT_tcam, IT_tlin =  ITfunction(source,scen,irt,2029)
        IT_2029 = ITrandom[case]*IT_tcam+(1.-ITrandom[case])*IT_tlin
      #
      #--store intensity reduction with peak constraint
      List_ChinaIRT[source][scen].append(irt) 
      #
      #--CO2 emissions when we do not consider the peak:
      emi['China']['co2'][2030][source][scen]['nopeak']=emi['China']['co2'][NDC['China']['baseline']]/            \
                                                            GDPrev[source][scen]['China'][NDC['China']['baseline']]* \
                                                            GDPrev[source][scen]['China'][2030]*(1.+NDC['China']['value'])
      #--CO2 emissions when we do consider the peak:
      emi['China']['co2'][2030][source][scen]['peak']=emi['China']['co2'][NDC['China']['baseline']]/              \
                                                            GDPrev[source][scen]['China'][NDC['China']['baseline']]* \
                                                            GDPrev[source][scen]['China'][2030]*(1.+irt)
      #
      #--third, in both cases, the CO2eq/CO2 ratio still has to be applied in order to get GHG emissions.
      for c in ['nopeak','peak']:
        #--Calculating China ratio CO2eq/CO2 for 2030
        ratioChina=ratioCfit[0][0]*np.log(emi['China']['co2'][2030][source][scen][c])+ratioCfit[0][1]+emiChinanonco2
        #--Applying ratio for China emissions in 2030.
        emi['China']['GHGexclLU'][2030][source][scen][c]=emi['China']['co2'][2030][source][scen][c]*ratioChina
        #
        #--fourth, we add LULUCF emissions:
        #--The Chinese NDC indicates a target of increasing the forest stock volume by 4.5 billion cubic meters in 2030
        #--compared to 2005. It also mentions an increase in 2014 compared to 2005 of 2.188 billion cubic meters.
        #--Bellassen et al. (2010) indicate that 1 cubic meter of forest is equivalent to removing 0.9175 tCO2.
        #--We therefore divide what is left of the 2030 target over the remaining fifteen years, and muliply it by 0.9175.
        #--The result is -139000 ktCO2, which we add to the 2010 CO2 LULUCF emissions.
        #
        emi['China']['co2lu'][2030][source][scen][c]=emi['China']['co2lu'][2010]-(4.5-2.188)*1000000.*0.9175/15.
        emi['China']['GHGtot'][2030][source][scen][c]=emi['China']['GHGexclLU'][2030][source][scen][c]\
                                                     +emi['China']['co2lu'][2030][source][scen][c]
      #
      #--Dealing with Other Annex I countries, Other Emerging countries and countries with NDC from Rest of World
      #
      for country in OtherAnnex1+OtherEmerging+RestOfWorldNDC:
        if NDC[country]['type']=='absolute':
          emi[country]['GHGtot'][2030][source][scen]=emi[country]['GHGtot'][NDC[country]['baseline']]*(1.+NDC[country]['value'])
        elif NDC[country]['type']=='BAU':
          emi[country]['GHGtot'][2030][source][scen]=NDC[country]['baseline']*(1.+NDC[country]['value'])
        elif NDC[country]['type']=='value':
          emi[country]['GHGtot'][2030][source][scen]=NDC[country]['value']
        else: #--intensity NDC types
          if country<>'Chile': #--For Singapore, the Philippines and Tunisia the target is applied to emissions including LULUCF following their NDCs
            emi[country]['GHGtot'][2030][source][scen]=emi[country]['GHGtot'][NDC[country]['baseline']]/  \
                                                       GDPrev[source][scen][country][NDC[country]['baseline']]*        \
                                                       GDPrev[source][scen][country][2030]*(1.+NDC[country]['value'])
          else:
            emi[country]['GHGexclLU'][2030][source][scen]=emi[country]['GHGexclLU'][NDC[country]['baseline']]/  \
                                                       GDPrev[source][scen][country][NDC[country]['baseline']]*        \
                                                       GDPrev[source][scen][country][2030]*(1.+NDC[country]['value'])
            #--Chile NDC has a separate target for LULUCF: 600+900+1200 ktCO2 in annual sequestration and reduction by 2030
            #--We add this number to the 2010 LULUCF emissions
            emi[country]['co2lu'][2030][source][scen]=emi[country]['co2lu'][2010]-2700
            emi[country]['GHGtot'][2030][source][scen]=emi[country]['GHGexclLU'][2030][source][scen]\
                                                      +emi[country]['co2lu'][2030][source][scen]
      #
      emi['OtherAnnex1']['GHGtot'][2030][source][scen]=np.sum(emi[country]['GHGtot'][2030][source][scen] for country in OtherAnnex1)
      emi['OtherEmerging']['GHGtot'][2030][source][scen]=np.sum(emi[country]['GHGtot'][2030][source][scen] for country in OtherEmerging)
      emi['RestOfWorldNDC']['GHGtot'][2030][source][scen]=np.sum(emi[country]['GHGtot'][2030][source][scen] for country in RestOfWorldNDC)
      #
      #--Dealing with Other high income oil exporters and the rest of Rest of World
      #
      #--Calculating 2030 emissions without LULUCF for these groups of countries, based on our hypothesis on intensity targets.
      for group in ['OtherOil','RestOfWorldNoNDC']:
        for country in eval(group):
          if GDPrev[source][scen][country][2005]<>0.:
            emi[country]['GHGexclLU'][2030][source][scen]=emi[country]['GHGexclLU'][NDC[group]['baseline']]/   \
                                                          GDPrev[source][scen][country][NDC[group]['baseline']]*  \
                                                          GDPrev[source][scen][country][2030]*(1.+NDC[group]['value'])
          else:#--when GDP values are missing, we use the ratio of the sums of GDP for countries in the group between baseline and target years.
            emi[country]['GHGexclLU'][2030][source][scen]=emi[country]['GHGexclLU'][NDC[group]['baseline']]/   \
                                                          np.sum(GDPrev[source][scen][country][NDC[group]['baseline']] for country in eval(group))*   \
                                                          np.sum(GDPrev[source][scen][country][2030] for country in eval(group))*(1.+NDC[group]['value'])
      #
      #--For countries with intensity targets, LULUCF emissions still have to be added afterwards; 
      #--the intensity target is applied to GHG emissions without LULUCF.
      #--Hypothesis: LULUCF emissions in 2030 at 2010 level. 
      for country in eval('OtherOil'):
        emi[country]['co2lu'][2030][source][scen]=emi[country]['co2lu'][2010]
      #
      #--Reduction in LULUCF emissions of RestOfWorld in 2030 compared to 2010
      for country in eval('RestOfWorldNoNDC'):
        emi[country]['co2lu'][2030][source][scen]=emi[country]['co2lu'][2010]*(1.-ROWredLU)
      #
      #--adding GHGexclLU + co2lu
      for group in ['OtherOil','RestOfWorldNoNDC']:
        for country in eval(group):
          emi[country]['GHGtot'][2030][source][scen]=emi[country]['GHGexclLU'][2030][source][scen]+emi[country]['co2lu'][2030][source][scen]
        emi[group]['GHGtot'][2030][source][scen]=np.sum(emi[country]['GHGtot'][2030][source][scen] for country in eval(group))
      emi['RestOfWorld']['GHGtot'][2030][source][scen]=emi['RestOfWorldNDC']['GHGtot'][2030][source][scen]+emi['RestOfWorldNoNDC']['GHGtot'][2030][source][scen]
      #
      #--Preparing output
      for country in WorldUnique+WorldOtherGroups:
        if country<>'China':
          List_Emi[source][scen][country]['GHGtot'][2030].append(emi[country]['GHGtot'][2030][source][scen])
        else:
          List_Emi[source][scen][country]['GHGtot'][2030].append(emi[country]['GHGtot'][2030][source][scen]['peak'])
      for emtype in ['co2']+typeOther:
        for c in ['nopeak','peak']:
          List_ChinaEmi[source][scen][emtype][2030][c].append(emi['China'][emtype][2030][source][scen][c])
      #
      #--Summing emissions for 2030
      WorldEmi2030=np.sum(emi[country]['GHGtot'][2030][source][scen] for country in MainCountries+Transport+WorldOtherGroups if country<>'China')
      WorldEmi2030+=emi['China']['GHGtot'][2030][source][scen]['peak']
      List_WorldEmi[source][scen]['GHGtot'][2030].append(WorldEmi2030/1.e6)
    #
#
#----------------------------------------------------------------------------------------------------------------------------------------
#--Deallocate some intermediate variables 
#
del WorldEmi2030, ROWredLU, emiChinanonco2
for country in emi.keys():
  if country!='WorldUnique':
    for emtype in emi[country].keys():
      if 2030 in emi[country][emtype].keys(): del emi[country][emtype][2030]
for emtype in emi['WorldUnique'].keys(): 
  if 2030 in emi['WorldUnique'][emtype].keys(): del emi['WorldUnique'][emtype][2030]
for country in NDC.keys():
  del NDC[country]['value']
#
#----------------------------------------------------------------------------------------------------------------------------------------
#--Outputs
#
print 'List of ',len(WorldUnique), 'entries used for computing world total emissions=', WorldUnique
#
print 'World GHG emissions for 2010=', "%.2f" % emi['WorldUnique']['GHGtot'][2010],'GtCO2eq'
#
for scen in SSPscenarios:
  for source in GDPsources:
    print 'World GHG emissions for 2030 ', source.ljust(6), scen.ljust(5), '=', \
                           "%.2f" % np.average(List_WorldEmi[source][scen]['GHGtot'][2030]),\
                           '+/-', "%.2f" % np.std(List_WorldEmi[source][scen]['GHGtot'][2030]),'GtCO2eq'
#
for scen in SSPscenarios:
  for source in GDPsources:
    print 'Chinese GHG emissions for 2030 ', source.ljust(6), scen.ljust(5), '=', \
               "%.2f" % (np.average(List_ChinaEmi[source][scen]['co2'][2030]['nopeak'])/1.e6), \
               '+/-', "%.2f" % (np.std(List_ChinaEmi[source][scen]['co2'][2030]['peak'])/1.e6), 'GtCO2'
#
#----------------------------------------------------------------------------------------------------------------------------------------
#--For each country in WorldUnique and groupings in WorldOtherGroups one can check the values of the 2030 emissions 
# List_Emi[source][scen][country]['GHGtot'][2030]
# np.average(List_Emi[source][scen][country]['GHGtot'][2030]) for computing the average of the distribution
# np.std(List_Emi[source][scen][country]['GHGtot'][2030]) for computing the standard deviation of the distribution
#
# For China, one can check the values of the 2030 emissions with and without the peak constraint
# List_ChinaEmi[source][scen]['GHGtot'][2030]['peak'] and List_ChinaEmi[source][scen]['GHGtot'][2030]['nopeak']
#
#--Examples
source='OECD'
scen='SSP2'
#
country='India'
print country, source, scen, 'GHG emissions=', \
               "%.2f" % (np.average(List_Emi[source][scen][country]['GHGtot'][2030])/1.e6),   \
               '+/-', "%.2f" % (np.std(List_Emi[source][scen][country]['GHGtot'][2030])/1.e6), 'GtCO2eq'
#
country='China'
print country, source, scen, 'CO2 emissions=', \
               "%.2f" % (np.average(List_ChinaEmi[source][scen]['co2'][2030]['nopeak'])/1.e6), \
               '+/-', "%.2f" % (np.std(List_ChinaEmi[source][scen]['co2'][2030]['peak'])/1.e6), 'GtCO2'
#
country='China'
percentile=10.
print country, source, scen, str(percentile)+'th percentile of CO2 emissions=', \
               "%.2f" % (np.percentile(List_ChinaEmi[source][scen]['co2'][2030]['peak'],percentile)/1.e6), 'GtCO2'
#
histo, bin_edges=np.histogram(List_ChinaIRT[source][scen],bins=40,range=[-1.,-0.6])
print country, source, scen, 'histogram of Intensity reduction values for China: ', histo 
#

print("Computation complete")

#--END