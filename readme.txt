This folder contains:

1) Datasets S1 and S2 and the code used for their analysis, see the file « sse.R ». These datasets have been considered in ''F. Wong, J. J. Collins, Evidence that coronavirus superspreading is fat-tailed. PNAS 117, 29416-29418 (2020).''

2) Database S3, see the file « traceDatSaved.Rdata », and the code used for its analysis, see « sse india.R ». This database has been considered in ''C. Kremer, A. Torneri, S. Boesmans, H. Meuwissen, S. Verdonschot, K. Vanden Driessche, C. L. Althaus, C. Faes, N. Hens, Quantifying superspreading for COVID-19 using Poisson mixture distributions. Sci. Rep. 11, 14107 (2021).''

3) Database S4, see the file « sse_korea_2020.txt », and the code used for its analysis, see « sse korea 2020.R ». This dataset has been considered in ''F. Wong, J. J. Collins, Evidence that coronavirus superspreading is fat-tailed. PNAS 117, 29416-29418 (2020).''

4) Database S5, see the file « sse_korea_2021_period_2.txt », and the code used for its analysis, see « sse korea 2021.R ». This database has been considered in ''S. Ryu, D. Kim, J.-S. Lim, S. T. Ali, B. J. Cowling, Serial interval and transmission dynamics during SARS-CoV-2 Delta variant predominance, South Korea. Emerg. Infect. Dis. 28, 407-410 (2022).''

5) Databases S6 and S7, see the files « colorado.txt » (Colorado data) and « clusters.Rdata » (all data except Colorado), and the codes used for their analysis, see « fit cluster XXX.R ». These databases have been compiled from several sources:
  Colorado: data originally taken from https://covid19.colorado.gov/covid19-outbreak-data, see the file "COVID-19 OB Weekly Report 06 02 2021.xlsx" ("resolved" sheet). All clusters with 2 or more cases are reported.
  Hong Kong: data originally taken from https://www.chp.gov.hk/files/pdf/local_situation_covid19_en.pdf, see the file "Hong-Kong-local-situation-covid19-en.pdf" ("Latest situation of cases of COVID-19 (as of 6 September 2021)"). Only clusters with 10 or more cases are reported.
  Kerala: data originally taken from https://covid19jagratha.kerala.nic.in/home/clusterList, available on Kaggle at https://www.kaggle.com/datasets/sebinjames/covid-clusters-in-kerala?select=KeralaClusters_August21_Update.csv (file "Kerala-Clusters\_August21\_Up").
  Oregon: data originally taken from https://www.oregon.gov/oha/covid19/Documents/DataReports/493Weekly-Outbreak-COVID-19-Report-2021-08-25-FINAL.pdfutm_medium=email&utm_494source=govdelivery, files "Weekly-Outbreak-COVID-19-Report-2021-08-04-FINAL", "Weekly-Outbreak-COVID-19-Report-2021-08-11-FINAL", "Weekly-Outbreak-COVID-19-Report-2021-08-18-FINAL", "Weekly-Outbreak-COVID-19-Report-2021-08-25-FINAL" and "Weekly-Outbreak-COVID-19-Report-2021-09-01-FINAL" (Care Facility, Senior Living Communities and Congregate Living Settings Reports and Workplace resolved outbreaks). All the outbreaks with 3 ore more confirmed cases or one or more deaths are reported.
  Australia, Brazil, California, Canada, China, France, Italy, New Jersey, Singapore, South Korea and United Kingdom: data taken from https://kmswinkels.medium.com/covid-19-superspreading-events-database-4c0a7aa2342b, file "SARS-CoV-2 Superspreading Events from Around the World" (column "Total Cases").

6) The codes used for robustness checks implemented on the cluster size data, see « robustness_XXX.R » (robustness due to missing data) and « robustness_XXX_bis.R » (robustness due to poor recording).

7) The code used to produce Figure 1, called « hillplot_burr_fig1.R ».
