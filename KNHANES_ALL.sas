%macro _hrqol_text(value);%sysfunc(kcvt(%sysfunc(inputc(&value,$hex%length(&value).)),utf-8,%sysfunc(getoption(encoding))))%mend;
%let survey_name=KNHANES;

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

data cln.cleaned(compress=binary);
 set source.KNHANES_Data(keep=year age sex edu marri_1 marri_2 DI1_dg DI1_lt DE1_dg DE1_lt
   BD1_11 BS1_1 BS3_1 HE_BMI wt_itvex kstrata psu
   BE3_11 BE3_21 BE3_31 BE3_71 BE3_75 BE3_81 BE3_85
   DI1_2 DE1_31 DE1_32 DE1_33 LQ_1EQL LQ_2EQL LQ_3EQL LQ_4EQL LQ_5EQL region allownc OR1 BP1 DM4_dg DM4_lt DM1_dg DM1_lt DI2_dg DI2_lt DI6_dg DI6_lt DI5_dg DI5_lt DI3_dg DI3_lt);
 source_row_id+1;
 length strata_year cluster_year stratum_year $40 psu_year $80;
 array raw_eq[5] LQ_1EQL LQ_2EQL LQ_3EQL LQ_4EQL LQ_5EQL;
 array bin_eq[5] eq11 eq21 eq31 eq41 eq51;
 eq_valid=1;
 do dim=1 to 5;
   bin_eq[dim]=.;
   if raw_eq[dim]=1 then bin_eq[dim]=0;
   else if raw_eq[dim] in(2,3) then bin_eq[dim]=1;
   else eq_valid=0;
 end;
 hypertension=.;diabetes=.;
 if DI1_dg=1 then hypertension=1;
 else if DI1_dg=0 then hypertension=2;

 else if DI1_dg=8 and DI1_lt=0 then hypertension=2;
 if DE1_dg=1 then diabetes=1;
 else if DE1_dg=0 then diabetes=2;
 else if DE1_dg=8 and DE1_lt=0 then diabetes=2;
 married=9;
 if marri_1=2 then married=2;
 if marri_2=1 then married=1;
 else if marri_2 in(2,4) then married=3;
 else if marri_2=3 then married=4;
 edu_cat=9;
 if edu in(1,2) then edu_cat=1;
 else if edu in(3,4) then edu_cat=edu;
 drink=9;
 if BD1_11=8 then drink=1;
 else if BD1_11=1 then drink=2;
 else if BD1_11 in(2,3,4,5,6) then drink=3;
 smoking=9;
 if BS1_1=3 then smoking=0;
 if BS3_1=3 then smoking=1;
 else if BS3_1 in(1,2) then smoking=2;
 gao=9;zhong=9;bu=9;
 if year<=2013 then do;
   if BE3_11=1 then gao=0;else if BE3_11 in(2,3,4,5,6,7,8) then gao=1;
   if BE3_21=1 then zhong=0;else if BE3_21 in(2,3,4,5,6,7,8) then zhong=1;
 end;
 else do;
   if BE3_71=1 or BE3_75=1 then gao=1;
   else if BE3_71=2 and BE3_75=2 then gao=0;
   if BE3_81=1 or BE3_85=1 then zhong=1;
   else if BE3_81=2 and BE3_85=2 then zhong=0;
 end;
 if BE3_31=1 then bu=0;else if BE3_31 in(2,3,4,5,6,7,8) then bu=1;
 BMI=.;if HE_BMI>0 then BMI=HE_BMI;
 weight2=wt_itvex/15.5;
 if year=2007 then weight2=weight2*.5;
 strata_year=catx('_',put(year,4.),strip(vvalue(kstrata)));
 cluster_year=catx('_',put(year,4.),strip(vvalue(psu)));
 diabetes_insulin=.;diabetes_medication=.;diabetes_non_drug=.;hypertension_medication=.;
 if DE1_31=1 then diabetes_insulin=1;else if DE1_31=0 then diabetes_insulin=2;
 if DE1_32=1 then diabetes_medication=1;else if DE1_32=0 then diabetes_medication=2;
 if DE1_33=1 then diabetes_non_drug=1;else if DE1_33=0 then diabetes_non_drug=2;
 if DI1_2 in(1,2,3,4) then hypertension_medication=1;
 else if DI1_2=5 then hypertension_medication=2;
 age_valid=(not missing(age) and age>=18);
 disease_valid=(hypertension in(1,2) and diabetes in(1,2));
 bmi_valid=(not missing(BMI) and BMI>0);
 design_valid=(not missing(kstrata) and not missing(psu) and weight2>0);
 sex_valid=(sex in(1,2));
 covariates_valid=(married in(1,2,3,4) and edu_cat in(1,3,4) and drink in(1,2,3)
  and smoking in(0,1,2) and gao in(0,1) and zhong in(0,1) and bu in(0,1));
 model_eligible=(age_valid and eq_valid and disease_valid and bmi_valid and design_valid and sex_valid and covariates_valid);
 
 josa_year=year;stratum_year=strata_year;psu_year=cluster_year;
 livelihood_history=.;if allownc in(10,20) then livelihood_history=allownc;
 oral_health=.;if OR1 in(1,2,3,4,5) then oral_health=OR1;
 perceived_stress=.;if BP1 in(1,2,3,4) then perceived_stress=BP1;
 array diagnosis_raw[6] DM4_dg DM1_dg DI2_dg DI6_dg DI5_dg DI3_dg;
 array lifetime_raw[6] DM4_lt DM1_lt DI2_lt DI6_lt DI5_lt DI3_lt;
 array disease_extra[6] osteoporosis arthritis dyslipidemia angina myocardial_infarction stroke;
 do dim=1 to 6;
  disease_extra[dim]=.;
  if diagnosis_raw[dim]=1 then disease_extra[dim]=1;
  else if diagnosis_raw[dim]=0 or (diagnosis_raw[dim]=8 and lifetime_raw[dim]=0) then disease_extra[dim]=2;
 end;
 sens1=model_eligible and nmiss(livelihood_history,oral_health,perceived_stress)=0;
 sens2=model_eligible and not missing(osteoporosis);
 sens3=model_eligible and nmiss(arthritis,dyslipidemia)=0;
 sens4=model_eligible and nmiss(angina,myocardial_infarction)=0;
 sens5=model_eligible and not missing(stroke);
 array old_regions[16] _temporary_ (11,26,27,28,29,30,31,41,42,43,44,45,46,47,48,50);
 array new_regions[17] _temporary_ (11,26,27,28,29,30,31,36,41,42,43,44,45,46,47,48,50);
 city_code_num=.;
 if year<=2015 and 1<=region and region<=16 then city_code_num=old_regions[region];
 else if year>2015 and 1<=region and region<=17 then city_code_num=new_regions[region];
 %region_names;
 source_record=1;
 valid_demographic=sex_valid and married in(1,2,3,4) and edu_cat in(1,3,4);
 valid_behavior=smoking in(0,1,2) and drink in(1,2,3);
 valid_activity=gao in(0,1) and zhong in(0,1) and bu in(0,1);
 drop dim;
run;

data cln.primary;set cln.cleaned;if model_eligible;run;
data cln.descriptive;
 set cln.cleaned;
 if age_valid and eq_valid and disease_valid and weight2>0;
run;

%flow_statistics(data=cln.cleaned,
 conditions=source_record age_valid eq_valid disease_valid bmi_valid valid_demographic valid_behavior valid_activity design_valid,
 keys=source|adult|eq|disease|bmi|demographic|behavior|activity|design,
 labels=Source|Adult age|EQ-5D responses|Disease responses|Anthropometry|Demographics|Smoking and alcohol|Physical activity|Survey design);
%baseline_statistics(data=cln.primary,
 categorical=age_group sex married edu_cat hypertension diabetes eq11 eq21 eq31 eq41 eq51 smoking drink gao zhong bu);
%regional_statistics(data=cln.descriptive);
data work.eq_city_year;set res.city_year;if indicator in('EQ1','EQ2','EQ3','EQ4','EQ5');run;
data work.eq_city_pooled;set res.city_pooled;if indicator in('EQ1','EQ2','EQ3','EQ4','EQ5');run;
data work.disease_city_year;set res.city_year;if indicator in('Hypertension','Diabetes');run;
data work.disease_city_pooled;set res.city_pooled;if indicator in('Hypertension','Diabetes');run;
%export_csv(data=work.eq_city_year,file=KNHANES_%_hrqol_text(E4BA94E7BBB4)_%_hrqol_text(E59CB0E58CBAE5B9B4E4BBBD));
%export_csv(data=work.eq_city_pooled,file=KNHANES_%_hrqol_text(E4BA94E7BBB4)_%_hrqol_text(E59CB0E58CBAE59088E5B9B6));
%export_csv(data=work.disease_city_year,file=KNHANES_%_hrqol_text(E9AB98E8A180E58E8BE7B396E5B0BFE79785)_%_hrqol_text(E59CB0E58CBAE5B9B4E4BBBD));
%export_csv(data=work.disease_city_pooled,file=KNHANES_%_hrqol_text(E9AB98E8A180E58E8BE7B396E5B0BFE79785)_%_hrqol_text(E59CB0E58CBAE59088E5B9B6));

%let base_covariates=married age sex edu_cat BMI smoking drink gao zhong bu;
%let base_classes=married(param=ref ref='1') edu_cat(param=ref ref='1')
 drink(param=ref ref='1') smoking(param=ref ref='1') gao(param=ref ref='1')
 zhong(param=ref ref='1') bu(param=ref ref='1') sex(param=ref ref='1');
%let survey_options=;
%init_results;
data work.main;set cln.primary;eligible=1;run;
%fit_model(data=work.main,analysis=main,scheme=main);
%macro sensitivity(s=,extra=,extra_class=);
 data work.sensitivity;set cln.primary;if sens&s;eligible=1;run;
 %fit_model(data=work.sensitivity,analysis=sensitivity,scheme=&s,role=BASE);
 %fit_model(data=work.sensitivity,analysis=sensitivity,scheme=&s,
  extra=&extra,extra_class=&extra_class);
%mend;
%sensitivity(s=1,extra=livelihood_history oral_health perceived_stress,
 extra_class=livelihood_history(param=ref ref='20') oral_health(param=ref ref='1') perceived_stress(param=ref ref='4'));
%sensitivity(s=2,extra=osteoporosis,extra_class=osteoporosis(param=ref ref='2'));
%sensitivity(s=3,extra=arthritis dyslipidemia,extra_class=arthritis(param=ref ref='2') dyslipidemia(param=ref ref='2'));
%sensitivity(s=4,extra=angina myocardial_infarction,extra_class=angina(param=ref ref='2') myocardial_infarction(param=ref ref='2'));
%sensitivity(s=5,extra=stroke,extra_class=stroke(param=ref ref='2'));

%macro sampling;
 %local i p;
 %do i=1 %to 3;
  %let p=%scan(90 70 50,&i);
  data work.fixed_ids;
   infile "&fixed_input\KNHANES_sample_ids_&p..csv" dsd dlm=',' firstobs=2 encoding='utf-8';
   input source_row_id year;keep source_row_id;
  run;
  proc sort data=work.fixed_ids;by source_row_id;run;
  data cln.sample_&p;
   merge cln.primary(in=a) work.fixed_ids(in=b);by source_row_id;
   if a and b;eligible=1;
  run;
  %fit_model(data=cln.sample_&p,analysis=sampling,scheme=&p);
 %end;
%mend;
%sampling;

%macro treatment;
 %local i factor diagnosis other;
 %do i=1 %to 4;
  %let factor=%scan(diabetes_insulin diabetes_medication diabetes_non_drug hypertension_medication,&i);
  %let diagnosis=diabetes;%let other=hypertension;
  %if &i=4 %then %do;%let diagnosis=hypertension;%let other=diabetes;%end;
  data work.treatment;
   set cln.primary;
   if &diagnosis=1 and &factor in(1,2);eligible=1;
  run;
  %fit_model(data=work.treatment,analysis=treatment,scheme=&factor,factor=&factor,other=&other);
 %end;
%mend;
%treatment;
data work.sensitivity_estimates;set res.estimates;if analysis='sensitivity';run;
data work.sensitivity_groups;
 set res.groups;if analysis='sensitivity' and model_role='ADJUSTED';
run;
%export_csv(data=work.sensitivity_estimates,file=KNHANES_%_hrqol_text(E4BA94E7BB84E6958FE6849FE680A7)_%_hrqol_text(E4BCB0E8AEA1));
%export_csv(data=work.sensitivity_groups,file=KNHANES_%_hrqol_text(E4BA94E7BB84E6958FE6849FE680A7)_%_hrqol_text(E4BABAE695B0));
%finish_results(prefix=KNHANES);
