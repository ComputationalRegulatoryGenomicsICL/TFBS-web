I am sorry, but documentation for TFBS modules is STILL not available. 
However, the modules themselves have evolved to 0.05 version.

New things in 0.05:

* pwmsearch C code now implemented as a proper XS extension -
  pwm_searchPFF.c and pwm_search.h are no longer used
* added pattern generator (TFBS::PatternGen::Gibbs) that works with 
  Chip Lawrence's Gibbs program
* added database creation and storage methods to JASPAR2 database
  interface (TFBS::DB::JASPAR2)
* added a number of methods to Matrix interfaces
* Fixed a number of 0.01 version bugs (and probably introduced some
  new ones in the process)

If you have any questions until the documantation is ready, please 
contact me at Boris.Lenhard@cgb.ki.se ;  I will be happy to help you.

The TFBS code can be distributed under the same terms as perl itself.

/Boris