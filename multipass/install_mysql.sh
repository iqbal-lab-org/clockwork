#!/usr/bin/env bash
DEBIAN_FRONTEND=noninteractive apt-get -y install mysql-server
mysql -e "CREATE USER 'test_user'@'localhost' IDENTIFIED BY 'test_password'"
mysql -e "GRANT ALL PRIVILEGES ON *.* TO 'test_user'@'localhost'"
mysql -e "FLUSH PRIVILEGES"
