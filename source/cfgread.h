//-------------------------------------------------------------------------------
//  Configuration Data Read
//-------------------------------------------------------------------------------

// default configuration data file
// if exist under the current directory, default variables in C-REX can be changed

static std::string CONFIG_FILE = "config.dat";

const int WORD_LENGTH = 256;
bool CFGRead(const std::string, char *);
