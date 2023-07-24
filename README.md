# GreenMO_Artefacts
Artefact submission for Mobicom'23 paper #121, "GreenMO"

## Data
The radio traces collected by us (iq-data) needs to be downloaded to run these evaluations. Please use ["FDT"](http://monalisa.cern.ch/FDT/download.html) to make the data transfer faster (usually takes about 5 minutes). You would need to setup Java JDK on your machine to run FDT (usually pre-installed on linux, for windows follow this [link](https://www.oracle.com/java/technologies/downloads/#jdk20-windows)). Open the command window and navigate to location where you have the "fdt.jar" file and run the following command to get the data required for artefact evaluation of GreenMO:


>  ```java -jar fdt.jar -p 54323 -pull -r -c 137.110.115.35 -r -d .\ag_data2\. /home/**/ag_data2/'''

Alternately, you can download from google drive: link. There might be some issues with local firewalls which may prevent binding to the FDT port.


## Trace Level Evaluations of GreenMO vs DBF/HBF (Fig. 9c, d)

Codes in reposit

>  
