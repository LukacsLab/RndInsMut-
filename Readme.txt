This code was developed by Julia Gala de Pablo
Some of the code was developed by Elton Rosas de Vasconcelo and edited by Julia Gala de Pablo
This code is very specific for our applications (to be published) of insertional mutagenesis, use at your own peril.
This code will be edited at publication date to make it user friendly for anyone that wants to reproduce our methods.

Example of our paired end input:
For a short insert before adapters ligation
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIISSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
ISSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
For a long insert before adapters ligation
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIISSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
ISSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
Where I is the adapter that needs removing (in our case it is part of the insertion), S is the sequence of interest and R is the reverse adapter. On R2 the adapter from the insertion is a single base, then the sequence follows and we may see the reverse of the adapter if the sequence is short.

That said, if you are working on something similar to our system and struggling, my inbox is open for discussion! 
