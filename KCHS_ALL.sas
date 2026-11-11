%macro _hrqol_text(value);%sysfunc(kcvt(%sysfunc(inputc(&value,$hex%length(&value).)),utf-8,%sysfunc(getoption(encoding))))%mend;
%let survey_name=KCHS;

%macro folders;
%global project_root source_root raw_input fixed_input output_parent run_name run_root;
%let project_root=D:\A %_hrqol_text(E68A95E7A8BF);
%if %sysfunc(sysexist(HRQOL_PROJECT_ROOT_HEX)) %then
 %let project_root=%_hrqol_text(%sysget(HRQOL_PROJECT_ROOT_HEX));
%else %if %sysfunc(sysexist(HRQOL_PROJECT_ROOT)) %then
 %let project_root=%sysget(HRQOL_PROJECT_ROOT);
%let source_root=D:\A %_hrqol_text(E68A95E7A8BF);
%if %sysfunc(sysexist(HRQOL_SOURCE_ROOT_HEX)) %then
 %let source_root=%_hrqol_text(%sysget(HRQOL_SOURCE_ROOT_HEX));
%else %if %sysfunc(sysexist(HRQOL_SOURCE_ROOT)) %then
 %let source_root=%sysget(HRQOL_SOURCE_ROOT);
%let raw_input=&source_root\%_hrqol_text(E58E9FE5A78BE695B0E68DAE);
%let fixed_input=&source_root\%_hrqol_text(E9A1B9E79BAEE8AEB0E5BD95E4B88EE6A0B8E69FA5)\%_hrqol_text(E695B0E68DAE)\%_hrqol_text(E59BBAE5AE9AE99A8FE69CBAE6A0B7E69CAC);
%let output_parent=&project_root\%_hrqol_text(E9A1B9E79BAEE8AEB0E5BD95E4B88EE6A0B8E69FA5)\%_hrqol_text(E7BB9FE8AEA1E8BF90E8A18CE7BB93E69E9C);
%let run_name=&survey_name._%sysfunc(datetime(),hex16.);
%let run_root=&output_parent\&run_name;
%if %sysfunc(sysexist(HRQOL_OUTPUT_ROOT_HEX)) %then
 %let run_root=%_hrqol_text(%sysget(HRQOL_OUTPUT_ROOT_HEX));
%else %if %sysfunc(sysexist(HRQOL_OUTPUT_ROOT)) %then
 %let run_root=%sysget(HRQOL_OUTPUT_ROOT);
%mend;
%folders;
options dlcreatedir validvarname=v7 compress=binary;
libname parent "&output_parent";
libname out "&run_root";
libname cln "&run_root\%_hrqol_text(E6B885E6B497E695B0E68DAE)";
libname res "&run_root\%_hrqol_text(E7BB9FE8AEA1E7BB93E69E9C)";
libname source "&raw_input" access=readonly;
libname fixed "&fixed_input" access=readonly;
data _null_;
 length path $1024;
 path=dcreate("%_hrqol_text(E7BB98E59BBEE695B0E68DAE)","&run_root");
 path=dcreate("%_hrqol_text(E68F8FE8BFB0E7BB9FE8AEA1)","&run_root");
run;
options nodlcreatedir;
filename analyslg "&run_root\%_hrqol_text(E58886E69E90E8BF90E8A18C).log" encoding='utf-8';
proc printto log=analyslg new;run;
ods html path="&run_root" file="SAS%_hrqol_text(E7BB9FE8AEA1E7BB93E69E9C).html" style=htmlblue;

%macro export_csv(data=,file=,folder=%_hrqol_text(E7BB98E59BBEE695B0E68DAE));
 data work.csv_output;
  set &data;
  format _numeric_ best32.;
  format _character_;
 run;
 filename csvout "&run_root\&folder\&file..csv" encoding='utf-8';
 proc export data=work.csv_output outfile=csvout dbms=csv replace;run;
 filename csvout clear;
%mend;

%macro init_results;
 data res.estimates;
  length analysis $16 scheme $64 model_role $16 exposure $48 used_years $180;
  length eq_dimension n_used beta se df estimate lower upper p 8;
  call missing(of _all_);
  stop;
 run;
 data res.groups;
  length analysis $16 scheme $64 model_role $16 exposure $48 group $3;
  length eq_dimension group_code n_group n_problem sum_weight weighted_event
         weighted_pct unweighted_pct 8;
  call missing(of _all_);
  stop;
 run;
%mend;

%macro fit_model(data=,analysis=,scheme=,role=ADJUSTED,extra=,extra_class=,factor=,other=);
 %local d rhs cls targets used_n used_years e target;
 %if %length(&factor) %then %do;
  %let rhs=&factor &other;
  %let cls=&factor(param=ref ref='2') &other(param=ref ref='2');
  %let targets=&factor;
 %end;
 %else %do;
  %let rhs=hypertension diabetes;
  %let cls=hypertension(param=ref ref='2') diabetes(param=ref ref='2');
  %let targets=hypertension diabetes;
 %end;
 proc sql noprint;
  select count(*) into :used_n trimmed from &data where eligible=1;
  select distinct josa_year into :used_years separated by ','
   from &data where eligible=1 order by josa_year;
 quit;
 %do d=1 %to 5;
  title "&survey_name &analysis &scheme &role EQ&d";
  ods output ParameterEstimates=work.pe OddsRatios=work.od;
  proc surveylogistic data=&data &survey_options namelen=64 varmethod=taylor;
   strata stratum_year;
   cluster psu_year;
   weight weight2;
   class &cls &base_classes &extra_class;
   model eq&d.1(event='1')=&rhs &base_covariates &extra;
  run;
  ods output close;
  data work.effects;
   length analysis $16 scheme $64 model_role $16 exposure $48 used_years $180;
   set work.pe;
   if strip(ClassVal0)='1';
   %if %length(&factor) %then %do;
    if strip(Variable)="&factor";
   %end;
   %else %do;
    if strip(Variable) in('hypertension','diabetes');
   %end;
   analysis="&analysis";scheme="&scheme";model_role="&role";
   exposure=strip(Variable);eq_dimension=&d;n_used=&used_n;
   beta=Estimate;se=StdErr;
   estimate=exp(beta);
   lower=exp(beta-tinv(.975,DF)*se);
   upper=exp(beta+tinv(.975,DF)*se);
   p=2*probt(-abs(beta/se),DF);
   used_years="&used_years";
   keep analysis scheme model_role exposure eq_dimension n_used beta se df
        estimate lower upper p used_years;
  run;
  proc append base=res.estimates data=work.effects force;run;
  %do e=1 %to %sysfunc(countw(&targets));
   %let target=%scan(&targets,&e);
   proc sql;
    create table work.group_counts as
     select &target as group_code,count(*) as n_group,
      sum(eq&d.1) as n_problem,sum(weight2) as sum_weight,
      sum(weight2*eq&d.1) as weighted_event
     from &data where eligible=1 group by &target;
   quit;
   data work.group_counts;
    length analysis $16 scheme $64 model_role $16 exposure $48 group $3;
    set work.group_counts;
    analysis="&analysis";scheme="&scheme";model_role="&role";
    exposure="&target";eq_dimension=&d;
    if group_code=1 then group='Yes';else if group_code=2 then group='No';
    weighted_pct=100*weighted_event/sum_weight;
    unweighted_pct=100*n_problem/n_group;
   run;
   proc append base=res.groups data=work.group_counts force;run;
  %end;
 %end;
%mend;

%macro region_names;
 length city_name_en $30;
 select(city_code_num);
 when(11) city_name_en='Seoul';when(26) city_name_en='Busan';
 when(27) city_name_en='Daegu';when(28) city_name_en='Incheon';
 when(29) city_name_en='Gwangju';when(30) city_name_en='Daejeon';
 when(31) city_name_en='Ulsan';when(36) city_name_en='Sejong';
 when(41) city_name_en='Gyeonggi-do';when(42) city_name_en='Gangwon-do';
 when(43) city_name_en='Chungcheongbuk-do';when(44) city_name_en='Chungcheongnam-do';
 when(45) city_name_en='Jeollabuk-do';when(46) city_name_en='Jeollanam-do';
 when(47) city_name_en='Gyeongsangbuk-do';when(48) city_name_en='Gyeongsangnam-do';
 when(50) city_name_en='Jeju-do';otherwise city_name_en='';end;
%mend;

%macro regional_statistics(data=);
 data work.region_long/view=work.region_long;
  set &data;
  length indicator $16 dimension $24 database $8 branch $16;
  array eq[5] eq11 eq21 eq31 eq41 eq51;
  array names[5] $24 _temporary_
   ('Mobility','Self-care','Usual activities','Pain/discomfort','Anxiety/depression');
  if not missing(city_name_en);
  if weight2<=0 then weight2=.;
  database="&survey_name";branch='descriptive';year=josa_year;
  do d=1 to 7;
   if d<=5 then do;indicator=cats('EQ',d);dimension=names[d];event=eq[d];end;
   else if d=6 then do;indicator='Hypertension';dimension='Hypertension';event=(hypertension=1);end;
   else do;indicator='Diabetes';dimension='Diabetes';event=(diabetes=1);end;
   if not missing(event) then output;
  end;
  keep database branch year josa_year city_code_num city_name_en indicator dimension event weight2;
 run;
 data work.region_values/view=work.region_values;
  set work.region_long;
  weighted_value=weight2*event;sq_weight=weight2*weight2;
 run;
 proc summary data=work.region_values nway;
  class database branch year josa_year city_code_num city_name_en indicator dimension;
  var event weight2 sq_weight weighted_value;
  output out=res.city_year(drop=_type_ _freq_)
   n(event)=n_valid sum(event weight2 sq_weight weighted_value)=n_event sum_weight sum_weight_sq weighted_event;
 run;
 data res.city_year;
  set res.city_year;
  weighted_pct=100*weighted_event/sum_weight;
  unweighted_pct=100*n_event/n_valid;
 run;
 %if &survey_name=KCHS %then %do;
  proc sql;
   create table work.region_grid as
    select y.database,y.branch,y.year,y.josa_year,c.city_code_num,c.city_name_en,
     v.indicator,v.dimension
    from (select distinct database,branch,year,josa_year from res.city_year) as y,
         (select distinct city_code_num,city_name_en from res.city_year) as c,
         (select distinct indicator,dimension from res.city_year) as v;
   create table work.complete_city_year as
    select g.*,coalesce(a.n_valid,0) as n_valid,coalesce(a.n_event,0) as n_event,
     coalesce(a.sum_weight,0) as sum_weight,coalesce(a.sum_weight_sq,0) as sum_weight_sq,
     coalesce(a.weighted_event,0) as weighted_event,a.weighted_pct,a.unweighted_pct
    from work.region_grid as g left join res.city_year as a
    on g.year=a.year and g.city_code_num=a.city_code_num and g.indicator=a.indicator;
  quit;
  data res.city_year;set work.complete_city_year;run;
 %end;
 proc sql;
  create table res.city_pooled as
   select database,branch,city_code_num,city_name_en,indicator,dimension,
    sum(n_valid) as n_valid,sum(n_event) as n_event,sum(sum_weight) as sum_weight,
    sum(weighted_event) as weighted_event,
    100*calculated weighted_event/calculated sum_weight as weighted_pct,
    100*calculated n_event/calculated n_valid as unweighted_pct
   from res.city_year group by database,branch,city_code_num,city_name_en,indicator,dimension;
  create table res.national_year as
   select database,branch,year,josa_year,indicator,dimension,
    sum(n_valid) as n_valid,sum(n_event) as n_event,sum(sum_weight) as sum_weight,
    sum(weighted_event) as weighted_event,
    100*calculated weighted_event/calculated sum_weight as weighted_pct,
    100*calculated n_event/calculated n_valid as unweighted_pct
   from res.city_year group by database,branch,year,josa_year,indicator,dimension;
 quit;
 %if &survey_name=KCHS %then %do;
  data work.city_export;set res.city_year;drop year;run;
  %export_csv(data=work.city_export,file=%_hrqol_text(E59F8EE5B882E5B9B4E4BBBD));
 %end;
 %else %do;
  %export_csv(data=res.city_year,file=%_hrqol_text(E59F8EE5B882E5B9B4E4BBBD));
 %end;
 %export_csv(data=res.city_pooled,file=%_hrqol_text(E59F8EE5B882E59088E5B9B6));
 %export_csv(data=res.national_year,file=%_hrqol_text(E5B9B4E4BBBDE695B0E68DAE)&survey_name,folder=%_hrqol_text(E68F8FE8BFB0E7BB9FE8AEA1));
%mend;

%macro baseline_statistics(data=,categorical=);
 %local expanded;
 %let expanded=&categorical current_smoking current_drinking vigorous_active moderate_active walking_active;
 data work.baseline_groups/view=work.baseline_groups;
  length group $8;
  set &data;
  age_group=1+(age>=30)+(age>=40)+(age>=50)+(age>=60)+(age>=70);
  %if &survey_name=KCHS %then %do;
   current_smoking=(smoke2=2);current_drinking=current_drinking_12m;
   vigorous_active=.;moderate_active=.;walking_active=.;
   if exercise1 in(0,3,7) then vigorous_active=(exercise1 in(3,7));
   if exercise2 in(0,3,7) then moderate_active=(exercise2 in(3,7));
   if exercise3 in(0,3,7) then walking_active=(exercise3 in(3,7));
  %end;
  %else %do;
   current_smoking=(smoking=2);current_drinking=(drink=3);
   vigorous_active=gao;moderate_active=zhong;walking_active=bu;
  %end;
  group='overall';output;
  if hypertension=1 then group='ht_yes';else group='ht_no';output;
  if diabetes=1 then group='dm_yes';else group='dm_no';output;
 run;
 %local i variable;
 proc summary data=work.baseline_groups missing;
  class group &expanded;
  types
  %do i=1 %to %sysfunc(countw(&expanded));
   group*%scan(&expanded,&i)
  %end;
  ;
  var weight2;
  output out=work.category_sums(drop=_type_ _freq_) n=n sum=weighted_sum;
 run;
 data work.category_sums;
  set work.category_sums;
  length variable $40;
  array values[*] &expanded;
  do i=1 to dim(values);
   if not missing(values[i]) then do;
    variable=vname(values[i]);category=values[i];output;
   end;
  end;
  keep group variable category n weighted_sum;
 run;
 proc sql;
  create table work.baseline_totals as
   select group,count(*) as denominator_n,sum(weight2) as denominator_weight
    from work.baseline_groups group by group;
  create table work.valid_totals as
   select group,variable,sum(n) as valid_n,sum(weighted_sum) as valid_weight
    from work.category_sums group by group,variable;
  create table res.baseline_categories as
   select a.variable,a.category,a.group,a.n,a.weighted_sum,
    case when a.variable in('vigorous_active','moderate_active','walking_active') then c.valid_n
     else b.denominator_n end as denominator_n,
    case when a.variable in('vigorous_active','moderate_active','walking_active') then c.valid_weight
     else b.denominator_weight end as denominator_weight,
    100*a.weighted_sum/calculated denominator_weight as weighted_pct,
    100*a.n/calculated denominator_n as unweighted_pct
   from work.category_sums as a left join work.baseline_totals as b on a.group=b.group
    left join work.valid_totals as c on a.group=c.group and a.variable=c.variable;
  create table res.baseline_continuous as
   select group,'age' as variable length=40,count(*) as n,
    sum(weight2) as denominator_weight,
    sum(weight2*age)/sum(weight2) as weighted_mean,
    sqrt(max(0,sum(weight2*age*age)/sum(weight2)-calculated weighted_mean**2)) as weighted_sd
   from work.baseline_groups group by group
   union all
   select group,'BMI' as variable,count(*) as n,sum(weight2) as denominator_weight,
    sum(weight2*BMI)/sum(weight2) as weighted_mean,
    sqrt(max(0,sum(weight2*BMI*BMI)/sum(weight2)-(sum(weight2*BMI)/sum(weight2))**2)) as weighted_sd
   from work.baseline_groups group by group;
 quit;
 %export_csv(data=res.baseline_categories,file=&survey_name._%_hrqol_text(E59FBAE7BABFE58886E7B1BBE7BB9FE8AEA1),folder=%_hrqol_text(E68F8FE8BFB0E7BB9FE8AEA1));
 %export_csv(data=res.baseline_continuous,file=&survey_name._%_hrqol_text(E59FBAE7BABFE8BF9EE7BBADE58F98E9878F),folder=%_hrqol_text(E68F8FE8BFB0E7BB9FE8AEA1));
 data work.table_cells;
  length category_text $24 display $60;
  set res.baseline_categories(in=c) res.baseline_continuous;
  if c then do;
   category_text=strip(put(category,best12.));
   display=catx(' ',strip(put(n,comma14.)),cats('(',strip(put(weighted_pct,6.1)),'%)'));
  end;
  else do;
   category_text='continuous';
   display=catx(' ',strip(put(weighted_mean,8.1)),'+/-',strip(put(weighted_sd,8.1)));
  end;
 run;
 proc sort data=work.table_cells;by variable category_text;run;
 proc transpose data=work.table_cells out=res.table1_display(drop=_name_);
  by variable category_text;id group;var display;
 run;
 %export_csv(data=res.table1_display,file=&survey_name._Table1_%_hrqol_text(E68E92E78988E8BE93E585A5),folder=%_hrqol_text(E68F8FE8BFB0E7BB9FE8AEA1));
%mend;

%macro flow_statistics(data=,conditions=,keys=,labels=);
 data work.flow_long/view=work.flow_long;
  set &data;
  length stage $40 label_en $60;
  array ok[*] &conditions;
  retained=1;
  do step=1 to dim(ok);
   n_enter=retained;
   retained=retained and ok[step];
   n_remaining=retained;n_excluded=n_enter-n_remaining;
   stage=scan("&keys",step,'|');label_en=scan("&labels",step,'|');output;
  end;
  keep josa_year step stage label_en n_enter n_excluded n_remaining;
 run;
 proc summary data=work.flow_long nway;
  class josa_year step stage label_en;var n_enter n_excluded n_remaining;
  output out=res.flow_annual(drop=_type_ _freq_) sum=;
 run;
 proc summary data=work.flow_long nway;
  class step stage label_en;var n_enter n_excluded n_remaining;
  output out=res.flow_total(drop=_type_ _freq_) sum=;
 run;
 %export_csv(data=res.flow_total,file=&survey_name._%_hrqol_text(E7BAB3E68E92E4BABAE695B0),folder=%_hrqol_text(E68F8FE8BFB0E7BB9FE8AEA1));
 %export_csv(data=res.flow_annual,file=&survey_name._%_hrqol_text(E7BAB3E68E92E98090E5B9B4E4BABAE695B0),folder=%_hrqol_text(E68F8FE8BFB0E7BB9FE8AEA1));
%mend;

%macro finish_results(prefix=);
 %export_csv(data=res.estimates,file=&prefix._%_hrqol_text(E4BCB0E8AEA1));
 %export_csv(data=res.groups,file=&prefix._%_hrqol_text(E4BABAE695B0));
 title "&survey_name model estimates";
 proc print data=res.estimates noobs;
  var analysis scheme model_role exposure eq_dimension n_used estimate lower upper p;
  format estimate lower upper 8.3 p pvalue8.4;
 run;
 title;
 ods html close;
 proc printto;run;
%mend;

%macro annual_design(year=,offset=,psu=,health=);
 proc import datafile="&raw_input\KCHS%_hrqol_text(E5B9B4E5BAA6)SPSS\KCHS_&year..sav"
  out=work.annual(keep=&psu &health) dbms=sav replace;run;
 data work.design_&year;
  length chs_psu $32 pbhlth_code $12;
  set work.annual(rename=(&psu=raw_psu &health=raw_health));
  source_row_id=_n_+&offset;
  chs_psu=strip(raw_psu);pbhlth_code=strip(raw_health);
  keep source_row_id chs_psu pbhlth_code;
 run;
%mend;

%macro design_2018(offset=);
filename sav2018 "&raw_input\KCHS%_hrqol_text(E5B9B4E5BAA6)SPSS\KCHS_2018.sav" recfm=n;
data work.design_2018;
 length packet cell $8 chs_psu $32 pbhlth_code $12;
 array codes[8] _temporary_;
 infile sav2018 recfm=n lrecl=32767;
 do k=1 to 56652;input byte $char1.;end;
 code_index=9;
 do case_index=1 to 228340;
  chs_psu='';pbhlth_code='';
  do storage_word=1 to 1112;
   code=0;
   do while(code=0);
    if code_index>8 then do;
     input packet $char8.;
     do k=1 to 8;codes[k]=rank(substr(packet,k,1));end;
     code_index=1;
    end;
    code=codes[code_index];code_index+1;
   end;
   cell='';
   if code=253 then input cell $char8.;
   if storage_word=13 then pbhlth_code=strip(cell);
   if storage_word=15 then chs_psu=strip(cell);
  end;
  source_row_id=&offset+case_index;output;
 end;
 stop;
 keep source_row_id chs_psu pbhlth_code;
run;

filename sav2018 clear;
%mend;
%annual_design(year=2011,offset=0,psu=jijum_cd,health=bogun_cd);
%annual_design(year=2012,offset=229226,psu=jijum_cd,health=bogun_cd);
%annual_design(year=2013,offset=458147,psu=jijum_cd,health=bogun_cd);
%annual_design(year=2014,offset=686928,psu=jijum_cd,health=bogun_cd);
%annual_design(year=2015,offset=915640,psu=jijum_cd,health=bogun_cd);
%annual_design(year=2016,offset=1144198,psu=jijum_cd,health=bogun_cd);
%annual_design(year=2017,offset=1372650,psu=jijum_cd,health=bogun_cd);
%design_2018(offset=1601031);
%annual_design(year=2019,offset=1829371,psu=spot_no,health=pbhlth_code);
%annual_design(year=2021,offset=2287739,psu=spot_no,health=pbhlth_code);
%annual_design(year=2022,offset=2516981,psu=spot_no,health=pbhlth_code);
%annual_design(year=2024,offset=2980518,psu=spot_no,health=pbhlth_code);
%annual_design(year=2025,offset=3212246,psu=spot_no,health=pbhlth_code);
data work.design;set work.design_2011 work.design_2012 work.design_2013 work.design_2014 work.design_2015 work.design_2016 work.design_2017 work.design_2018 work.design_2019 work.design_2021 work.design_2022 work.design_2024 work.design_2025;run;

data work.raw;
 set source.KCHS_Data;
 source_row_id=_n_;
run;
data cln.cleaned;
 merge work.raw(in=a) work.design;by source_row_id;if a;
 length stratum_year $40 psu_year $80;
 stratum_year=catx('_',put(josa_year,4.),strip(kstrata));
 psu_year=catx('_',put(josa_year,4.),pbhlth_code,chs_psu);
 raw_drink=drink2;
 valid_year=(josa_year in(2011,2012,2013,2014,2015,2016,2017,2018,2019,2021,2022,2024,2025));
 valid_age=(age>=19);
 valid_hw=(high>0 and weigh>0 and high not in(777,999,77777,88888,99999)
  and weigh not in(777,999,77777,88888,99999));
 descriptive_hw=(not missing(high) and not missing(weigh));
 BMI=.;if valid_hw then BMI=weigh/((high/100)**2);
 if mar not in(1,2,3,4,5,8) then mar=.;
 else if mar=2 then mar=4;else if mar=8 then mar=5;
 education=.;
 if edu in(1,2,3,4) then education=4;
 else if edu=5 then education=5;else if edu in(6,7,8) then education=8;
 if smoke2 not in(1,2,3,8) then smoke2=.;else if smoke2=1 then smoke2=2;
 if drink2 not in(1,2,3,4,5,6,8) then drink2=.;
 else if drink2 in(1,2,3,4) then drink2=5;else if drink2=6 then drink2=9;
 array exercise[3] exercise1 exercise2 exercise3;
 do i=1 to 3;
  if exercise[i] not in(0,1,2,3,4,5,6,7) then exercise[i]=.;
  else if exercise[i] in(1,2,3) then exercise[i]=3;
  else if exercise[i] in(4,5,6,7) then exercise[i]=7;
 end;
 if josa_year=2018 then do;exercise1=99;exercise2=99;end;
 if hypertension not in(1,2) then hypertension=.;
 if diabetes not in(1,2) then diabetes=.;
 if sex not in(1,2) then sex=.;
 array raw_eq[5] EQ1-EQ5;array bin_eq[5] eq11 eq21 eq31 eq41 eq51;
 valid_eq=1;
 do i=1 to 5;
  if raw_eq[i]=1 then bin_eq[i]=0;
  else if raw_eq[i] in(2,3) then bin_eq[i]=1;
  else do;bin_eq[i]=.;valid_eq=0;end;
 end;
 valid_disease=(nmiss(hypertension,diabetes)=0);
 valid_demographic=(nmiss(mar,education,age,sex)=0);
 valid_behavior=(nmiss(smoke2,drink2)=0);
 valid_activity=(nmiss(exercise1,exercise2,exercise3)=0);
 valid_design=(weight2>0 and not missing(kstrata) and not missing(chs_psu) and not missing(pbhlth_code));
 main_eligible=(valid_year and valid_age and valid_eq and valid_disease and valid_hw
  and valid_demographic and valid_behavior and valid_activity and valid_design);
 current_drinking_12m=.;
 if josa_year<=2019 then do;
  if raw_drink in(1,2,3,4,5) then current_drinking_12m=1;
  else if raw_drink=8 then current_drinking_12m=0;
 end;
 else do;
  if raw_drink in(2,3,4,5,6) then current_drinking_12m=1;
  else if raw_drink in(1,8) then current_drinking_12m=0;
 end;
 if fma_04z1 not in(1,2,3) then fma_04z1=.;
 if ora_01z1 not in(1,2,3,4,5) then ora_01z1=.;
 else if ora_01z1=2 then ora_01z1=1;else if ora_01z1=5 then ora_01z1=4;
 if mta_01z1 not in(1,2,3,4) then mta_01z1=.;else if mta_01z1=2 then mta_01z1=1;
 array disease_extra[6] osteoporosis arthritis dyslipidemia angina myocardial_infarction stroke;
 do i=1 to 6;if disease_extra[i] not in(1,2) then disease_extra[i]=.;end;
 array treatment[5] diabetes_insulin diabetes_medication diabetes_non_drug hypertension_medication hypertension_non_drug;
 do i=1 to 5;if treatment[i] not in(1,2) then treatment[i]=.;end;
 sens1=main_eligible and nmiss(fma_04z1,ora_01z1,mta_01z1)=0;
 sens2=main_eligible and josa_year in(2011,2012,2013,2016,2025) and not missing(osteoporosis);
 sens3=main_eligible and josa_year in(2011,2012,2013,2014,2015,2016,2017) and nmiss(arthritis,dyslipidemia)=0;
 sens4=main_eligible and josa_year in(2011,2012,2013) and nmiss(angina,myocardial_infarction)=0;
 sens5=main_eligible and josa_year in(2011,2012,2013,2014,2016) and not missing(stroke);
 city_code_num=input(substr(strip(city_id),1,2),best12.);
 if city_code_num=51 then city_code_num=42;if city_code_num=52 then city_code_num=45;
 %region_names;
 source_record=1;
 drop i;
run;
data cln.study_frame;set cln.cleaned;if valid_year;run;
data cln.primary;set cln.study_frame;if main_eligible;run;
data cln.descriptive;
 set cln.study_frame;
 if valid_eq and valid_disease and descriptive_hw;
run;

%flow_statistics(data=cln.cleaned,
 conditions=source_record valid_year valid_age valid_eq valid_disease valid_hw valid_demographic valid_behavior valid_activity valid_design,
 keys=source|years|adult|eq|disease|anthro|demographic|behavior|exercise|weight,
 labels=Source|Study years|Adult age|EQ-5D responses|Disease responses|Anthropometry|Demographics|Smoking and alcohol|Physical activity|Survey design);
%baseline_statistics(data=cln.primary,
 categorical=age_group sex mar education hypertension diabetes eq11 eq21 eq31 eq41 eq51 smoke2 drink2 current_drinking_12m exercise1 exercise2 exercise3);
%regional_statistics(data=cln.descriptive);

%let base_covariates=mar age sex education smoke2 drink2 exercise1 exercise2 exercise3 BMI;
%let base_classes=mar(param=ref ref='1') education(param=ref ref='4')
 exercise1(param=ref ref='3') exercise2(param=ref ref='3') exercise3(param=ref ref='3')
 smoke2(param=ref ref='2') drink2(param=ref ref='8');
%let survey_options=nomcar;
%init_results;
%macro kchs_model(analysis=,scheme=,role=ADJUSTED,years=,eligible=,extra=,extra_class=,factor=,other=);
 data work.model_frame;
  set cln.study_frame;
  if josa_year in(&years);
  eligible=(&eligible);
  if not eligible then call missing(eq11,eq21,eq31,eq41,eq51);
 run;
 %fit_model(data=work.model_frame,analysis=&analysis,scheme=&scheme,role=&role,
  extra=&extra,extra_class=&extra_class,factor=&factor,other=&other);
%mend;
%let years=2011 2012 2013 2014 2015 2016 2017 2018 2019 2021 2022 2024 2025;
%kchs_model(analysis=main,scheme=main,years=&years,eligible=main_eligible);
%macro sensitivity(s=,years=,extra=,extra_class=);
 %kchs_model(analysis=sensitivity,scheme=&s,role=BASE,years=&years,eligible=sens&s);
 %kchs_model(analysis=sensitivity,scheme=&s,years=&years,eligible=sens&s,extra=&extra,extra_class=&extra_class);
%mend;
%sensitivity(s=1,years=&years,extra=fma_04z1 ora_01z1 mta_01z1,
 extra_class=fma_04z1(param=ref ref='1') ora_01z1(param=ref ref='1') mta_01z1(param=ref ref='1'));
%sensitivity(s=2,years=2011 2012 2013 2016 2025,extra=osteoporosis,extra_class=osteoporosis(param=ref ref='1'));
%sensitivity(s=3,years=2011 2012 2013 2014 2015 2016 2017,extra=arthritis dyslipidemia,
 extra_class=arthritis(param=ref ref='1') dyslipidemia(param=ref ref='1'));
%sensitivity(s=4,years=2011 2012 2013,extra=angina myocardial_infarction,
 extra_class=angina(param=ref ref='1') myocardial_infarction(param=ref ref='1'));
%sensitivity(s=5,years=2011 2012 2013 2014 2016,extra=stroke,extra_class=stroke(param=ref ref='1'));

%macro sampling;
 %local p i;
 %do i=1 %to 3;
  %let p=%scan(90 70 50,&i);
  data work.fixed_ids;set fixed.KCHS_sample_ids_&p(keep=source_row_id);selected=1;run;
  proc sort data=work.fixed_ids;by source_row_id;run;
  data cln.sample_&p;
   merge cln.study_frame(in=a) work.fixed_ids;by source_row_id;
   if a;
   eligible=main_eligible and selected=1;
   if not eligible then call missing(eq11,eq21,eq31,eq41,eq51);
  run;
  %fit_model(data=cln.sample_&p,analysis=sampling,scheme=&p);
 %end;
%mend;
%sampling;

%macro treatment;
 %local i factor diagnosis other;
 %do i=1 %to 5;
  %let factor=%scan(diabetes_insulin diabetes_medication diabetes_non_drug hypertension_medication hypertension_non_drug,&i);
  %let diagnosis=diabetes;%let other=hypertension;
  %if &i>=4 %then %do;%let diagnosis=hypertension;%let other=diabetes;%end;
  %kchs_model(analysis=treatment,scheme=&factor,years=&years,
   eligible=main_eligible and &diagnosis=1 and &factor in(1,2),factor=&factor,other=&other);
 %end;
%mend;
%treatment;
%finish_results(prefix=CHS);
