********************************************************;
*
*	SAS code for simulation in
*	A Warning About Using Predicted Values to Estimate Descriptive Measures
*   Ross et al.
*	In press, AJE
*
*	Code by Rachael Ross, rkross@unc.edu
*
*******************************************************;

OPTIONS FORMCHAR="|----|+|---+=|-/\<>*"  symbolgen mprint mlogic notes
			source source2 msglevel=I mcompilenote=all mautosource;

*Generate data;
%macro gendata(sd,numsim=,sampsize=,outdat=);
data &outdat.;
call streaminit(13);
do j=1 to &numsim.;
	do i=1 to &sampsize.;
	sd = &sd.;

	*Observed GA at birth;
	Ey_z = rand('WEIBULL',31.61915,39.72913); 

	*True GA at birth;
	y = rand('normal',Ey_z, &sd./7);

	*Indicators for preterm and postterm;
	if Ey_z < 37 then obs_pretb = 1; else obs_pretb = 0;
	if y < 37 then true_pretb = 1; else true_pretb = 0;
	if Ey_z ge 42 then obs_posttb = 1; else obs_posttb = 0;
	if y ge 42 then true_posttb = 1; else true_posttb = 0;

	output;
end;
end;
run;
%mend;

%gendata(sd=2.7,numsim=5000,sampsize=2000,outdat=out1);
%gendata(sd=6.1,numsim=5000,sampsize=2000,outdat=out2);

*Analyze simulated data;
%macro analyze;
proc datasets nolist; 
		delete results ; 
	quit; run;

%let list=out1 out2;
%do t = 1 %to 2;
%let data = %scan(&list.,&t.);

proc sql noprint;
select mean(true_pretb), mean(true_posttb) into :true_pretb, :true_posttb from &data.; quit;

%let var=obs_pretb true_pretb obs_posttb true_posttb;
*Loop over vars;
%do s = 1 %to 4;
%let var_ = %scan(&var.,&s.);

%if &s. < 3 %then %do; %let true = &true_pretb.; %end;
%else %do; %let true = &true_posttb. ;%end;
%put &true.;

proc sql;
create table byj as 
select count(*) as n, 
	mean(&var_.) as mean, 
	sqrt(calculated mean*(1-calculated mean)/calculated n) as se,
	calculated mean - 1.96* calculated se as lcl,
	calculated mean + 1.96* calculated se as ucl
from &data. group by j; quit;

proc sql;
create table summary as
select "&data." as tag, "&var_." as var, &true. as truth,
	mean(mean) as mean, 
	calculated mean - &true. as bias, 
	std(mean) as ese, 
	mean(se) as amse,
	mean(case when lcl le &true. and ucl ge &true. then 1 else 0 end) as cover,
	count(*) as countj
from byj; quit;


%if %sysfunc(exist(results)) %then %do; 
	data results; length tag $10 var $15;
	set results summary;
	run;
	%end;
	%else %do;
	data results; set summary; run;
	%end;

%end;

%end;
%mend;

%analyze;
proc sort data=results;
by tag truth var;
run;
proc print data=results noobs;
run;


*Multiple imputation correction approach;
%macro mi(data,nimp);

*Get the truth;
proc sql noprint;
select mean(y), mean(true_pretb), mean(true_posttb) into :true_mean, :true_pretb, :true_posttb from &data.; quit;

*Multiply impute;
data imputed;
	call streaminit(15);
	set &data.;
	do t = 1 to &nimp.;

		delta = 1;
		draw = rand("normal",Ey_z,sd/7);
		if draw < 37 then new_pretb = 1; else new_pretb = 0;
		if draw ge 42 then new_posttb = 1; else new_posttb = 0;
	output;
	end;
	run;

*Analyze within imputations;
proc sql;
create table byimp0 as
select j, t, "new_mean" as Parm, &true_mean. as true, 
	mean(draw) as p, stderr(draw)**2 as var 
from imputed group by j, t; quit;

proc sql;
create table byimp1 as
select j, t, "new_pretb" as Parm, &true_pretb. as true, 
	mean(new_pretb) as p, calculated p*(1-calculated p)/count(*) as var 
from imputed group by j, t; quit;

proc sql;
create table byimp2 as
select j, t, "new_posttb" as Parm, &true_posttb. as true, 
	mean(new_posttb) as p, calculated p*(1-calculated p)/count(*) as var 
from imputed group by j, t; quit;


data byimp; 
length Parm $10.; 
set byimp0 byimp1 byimp2; 
run;

*Combine imputations;
proc sql;
create table byj as
select j, Parm, true, 
	mean(p) as estimate, 
	1/count(*)*sum(var) as within, 
	(1+1/count(*))*(1/(count(*)-1))*sum((var - meanvar)**2) as btw, 
	sqrt(calculated within + calculated btw) as se, 
	calculated estimate - 1.96*calculated se as lcl, calculated estimate + 1.96*calculated se as ucl 
from (select *, mean(var) as meanvar from byimp group by j, Parm) 
group by j, meanvar, Parm, true; 
quit;

*Summarize across replicates;
proc sql;
create table bymi_&data. as
select "&data." as tag, Parm as var, true as truth, 
	mean(estimate) as mean, 
	mean(estimate - true) as bias, 
	mean(se) as amse,
	mean(case when lcl le true and ucl ge true then 1 else 0 end) as cover, 
	count(*) as countj
from byj group by Parm, true; 
quit;
%mend;

%mi(out2,20);

proc print data=bymi_out2 noobs;
run;


