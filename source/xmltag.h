static inline void XMLTagOpen    (std::string);
static inline void XMLTagClose   (std::string);
static inline void XMLTagVal     (std::string, int);
static inline void XMLTagVal     (std::string, double);
static inline void XMLTagVal     (std::string, std::string);
static inline void XMLTagIndent  (int);

static int default_indent_level = 3;
static int indentlevel = 0;


/**********************************************************/
/*      Indent Level Modifier                             */
/**********************************************************/
void XMLSetDefaultIndentLevel(int n){
  default_indent_level = n;
}


/**********************************************************/
/*      XML TAG Output                                    */
/**********************************************************/
void XMLTagOpen(std::string tag)
{
  XMLTagIndent(0);
  cout << "<" << tag << ">" << endl;
  indentlevel++;
}

void XMLTagOpen(std::string tag, std::string attr)
{
  XMLTagIndent(0);
  cout << "<" << tag << " " << attr << ">" << endl;
  indentlevel++;
}

void XMLTagClose(std::string tag)
{
  indentlevel--;
  XMLTagIndent(0);
  cout << "</" << tag << ">" << endl;
  if(indentlevel<0) indentlevel = 0;
}

void XMLTagVal(std::string tag, int val)
{
  XMLTagIndent(0);
  cout << "<" << tag << "> " << val << " </" << tag << ">" << endl;
}

void XMLTagVal(std::string tag, double val)
{
  XMLTagIndent(0);
  cout << "<" << tag << "> " << val << " </" << tag << ">" << endl;
}

void XMLTagVal(std::string tag, std::string val)
{
  XMLTagIndent(0);
  cout << "<" << tag << "> " << val << " </" << tag << ">" << endl;
}

void XMLTagIndent(int n)
{
  for(int i=0 ; i<n + indentlevel ; i++){
    for(int j=0 ; j<default_indent_level ; j++) cout << " ";
  }
}
