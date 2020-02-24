antlrcpp::Any::~Any ()
{
  if (this->get_typeinfo() == typeid(std::vector<double>*)) {
    //    std::cout << " Deleting double vec\n";
    //std::vector<double> *vals = this->as<std::vector<double>*>();
    //for (size_t i = 1; i < vals->size (); i++)
    //  std::cout << (*vals)[i] << " ";
    //std::cout << std::endl;
    //delete vals;
    this->clear ();
  }
  else if (this->get_typeinfo() == typeid(std::vector<bool>*)) {
    //    std::cout << " ALLOCATED boo; vec\n";
    //std::vector<bool> *vals = this->as<std::vector<bool>*>();
    //delete vals;
    this->clear ();
  }
  else if (this->get_typeinfo() == typeid(std::vector<size_t>*)) {
    //    std::cout << " ALLOCATED boo; vec\n";
    //std::vector<size_t> *vals = this->as<std::vector<size_t>*>();
    //delete vals;
    this->clear ();
  }
  else if (this->get_typeinfo() == typeid(Matrix*)) {
    //    std::cout << " ALLOCATED Matrix\n";
    //Matrix *vals = this->as<Matrix*>();
    //delete vals;
    this->clear ();
  }
  else if (this->get_typeinfo() == typeid(std::string*)) {
    // fixme double free or corruption (out)
    //std::cout << " ALLOCATED string\n";
    //std::string *vals = this->as<std::string*>();
    //    delete vals;
    this->clear ();
  }
  else if (this->get_typeinfo() == typeid(nullptr) ||
	   this->get_typeinfo() == typeid(double) ||
	   this->get_typeinfo() == typeid(int) ||
	   //	   this->get_typeinfo() == typeid(std::char_traits<char>)
	   //	   this->get_typeinfo() == typeid(std::char_traits<char>)
	   //	   this->get_typeinfo() == typeid(basic_string)
	   //	   	   this->get_typeinfo() == typeid(basic_string)
	   strstr (get_typeinfo().name (), "basic_string")
	   ) {
    //std::cout << " ALLOCATED string\n";
  }
  else {
    std::cout << "NOT CAUGHT "
              << "(" << demanglep (this) << ") = \n\t\""
	      << this->get_typeinfo().name () << "\""
              << std::endl;
  }
}
