

Hello there 😄 Welcome to this Jupyter Book instance to document mainly a data analysis workflow for analysing data from the technique : **S**aturated **T**ransposition **A**nalysis in **Y**east. 

Bellow some guidelines to get started using this tool in your own project 👇

## To reproduce an example of  Jupyter book, and you are a Windows user:


```{toggle} First! Try to reproduce this! Click the buttom to reveal :) 

The following workflow should succeed using a miniconda powershell terminal on Windows 10:

`Open Anaconda Power Shell Promt`

0. `conda install -c eoas_ubc runjb`
1. `conda install git`
2. `git clone https://github.com/eoas-ubc/quantecon-mini-example.git`
3. `cd quantecon-mini-example`
4. `git checkout windows`
5. `conda env create -f environment_win.yml`
6. `conda activate wintest`
7. `cd mini_book`
8. `runjb docs`

After the build, view the html with (local deployment):

`start docs\_build\html\index.html`
```
- For your own jupyter book just respect the same folder structure of the example, but with your own content! 

## How to deploy your book online? 

- Add a folder in the root of the repo called : `.github`
  - Add a subfolder called : `workflows`
  - Create a file called : `build.yml` where the instructions to GitHub to deploy your book are given. The following `build.yml` is generic as long as your **deploy Key** is called : `DEPLOY_KEY`. 

```yaml
name: Test Build
on:
  push:
    branches-ignore:
      - master
jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout 🛎️
        uses: actions/checkout@v2
        with:
          persist-credentials: false

      - name: Setup Miniconda
        uses: goanpeca/setup-miniconda@v1
        with:
          auto-update-conda: true
          auto-activate-base: false
          miniconda-version: 'latest'
          python-version: 3.7
          environment-file: environment_win.yml
          activate-environment: wintest

      - name: Install jupyter_book
        shell: bash -l {0}
        run:  pip install jupyter-book --pre 


      - name: Build Jupyter Book
        shell: bash -l {0}
        run: jb build  mini_book/docs/

      - name: Install SSH Client 🔑
        uses: webfactory/ssh-agent@v0.2.0
        with:
          ssh-private-key: ${{ secrets.DEPLOY_KEY }}
        


      - name: Deploy 🚀
        uses: JamesIves/github-pages-deploy-action@releases/v3
        with:
          SSH: true
          BRANCH: gh-pages
          FOLDER: mini_book/docs/_build/html/
          
```
- Create a new branch called : 'gh-pages'
- A public and private key just dedicated to this Jupyter Book repo 
- The content of the private key should be copied to the `Secrets` tab in the Settings of the repo. 
    - Create an SSH key with sufficient access privileges. For security reasons, don't use your personal SSH key but set up a dedicated one for use in GitHub Actions. Instructions 👇
    - Steps:
     1. Open Git Bash.
     2. Insert in bash `ssh-keygen -t rsa -b 4096 -C "your_email@example.com"` This creates a new ssh key, using the provided email as a label.
     > Generating public/private rsa key pair.
    When you're prompted to "Enter a file in which to save the key," press Enter. This accepts the default file location.
    > Enter a file in which to save the key (/c/Users/you/.ssh/id_rsa):[Press enter]
    At the prompt, type a secure passphrase. For more information, see "Working with SSH key passphrases".
    > Enter passphrase (empty for no passphrase): [Press enter] 
    - Make sure you don't have a passphrase set on the private key.
    > Enter same passphrase again: [Press enter again]


- In your repository, go to the Settings > Secrets menu and create a new secret. In this example, we'll call it SSH_PRIVATE_KEY. Put the contents of the private SSH key file into the contents field.
- This key should start with `-----BEGIN ... PRIVATE KEY-----`, consist of many lines and ends with `-----END ... PRIVATE KEY-----`.
- SSH private key format
    - If the private key is not in the PEM format, you will see an Error loading key `"(stdin)": invalid format message`.

    - Use `ssh-keygen -p -f path/to/your/key -m pem` to convert your key file to PEM, but be sure to make a backup of the file first 

- Add an empty file `.nojekyll` to the root of the repo to prevent github on deploying a jekyll website.


