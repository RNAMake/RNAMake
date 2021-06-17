## Generating options file from a csv file

**1.Craete the options csv file** 

* Create a csv file called ``options.csv`` and include in the app's folder

The csv file columns should follow this format

name,group,struct,parameter,help,default,required,flag,type


column|type|use
-----|-----|-----
name|string|Name of the option
group|string|Name of the option group
struct|string|Name of the sub-struct that the option belongs to
parameter|string|Name of the parameter in the code (Usually same as the name)
help|string|Help message text for the option
default|string|Default value for parameter (make sure to include quotations for string type options)
required|boolean|Indicate if the option is required or not
flag|boolean|Indicate if the option is a flag
type|string|The type of the parameter



**2.Running the script**

* Run the following script

``python generate_options.py --app <APP_NAME>``

The script will create an options header file called ```<APP_NAME>_options.h```, make sure to include this file in the app ```.cpp\.cc``` file.

