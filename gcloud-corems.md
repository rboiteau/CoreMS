
## CoreMS on gcloud

### Installing and initializing gcloud CLI

1. Follow these instructions to install: [gcloud-cli-install-instructions](https://cloud.google.com/sdk/docs/install-sdk)

2. When prompted, log in with the Google account with access to the corems-gcloud project. 

3. Select the following project ID: corems-377723

4. For zone, select us-central1-c


### Starting a gcloud compute session on manatee

1. From the terminal, run the following command

    `gcloud compute ssh manatee --zone=us-central1-c`

2. If this is the first time you are accessing the gcloud machine, you will be prompted to create a passphrase for a new ssh access key.  
