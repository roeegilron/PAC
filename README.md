Attempt at faster PAC code  
==========================

Files: 
-------------

* `CallerRoutine()` - Main function - some changes here to pre compute values 
* `ModIndex_v3()` - Function that computes MI  
* `ExtractHGHFOOpenField,mat` - sample data 
* `eegfilt()` - filtering from eeglab   


Main changes used to speed up PAC computation: 
-------------

* Use logical indexing 
* Pre compute values for MI function 
* Largest dimension first arrays 
* Linearize `Comodulogram` variable for faster `parfor` performance 

Tried and failed / modest gains: 
-------------
* Vectorizing code by pre computing indexing arrays 

To do: 
-------------
* Try using `gpuArray` 
* Implement surragtes as option (stats) 



