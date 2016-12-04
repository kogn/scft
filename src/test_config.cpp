//main.cpp
#include "Config.h"
int main()
{
	int port;
	std::string ipAddress;
	std::string username;
	std::string password;
	const char ConfigFile[]= "cfg/config.txt"; 
	Config configSettings(ConfigFile);
	
	port = configSettings.Read("port",port);
	ipAddress = configSettings.Read("ipAddress", ipAddress);
	username = configSettings.Read("username", username);
	password = configSettings.Read("password", password);
	std::cout<<"port:"<<port<<std::endl;
	std::cout<<"ipAddress:"<<ipAddress<<std::endl;
	std::cout<<"username:"<<username<<std::endl;
	std::cout<<"password:"<<password<<std::endl;
	
	return 0;
}

