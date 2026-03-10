# Push this project to GitHub

The repo is initialized with an initial commit on branch `main`. To publish it on GitHub:

## 1. Create the repository on GitHub

1. Open **https://github.com/new**
2. Set **Repository name** to: `MetagenomeGenerator` (or another name you prefer)
3. Set **Description** (optional): `Download NCBI genomes and build simulated metagenome FASTAs for training classifiers`
4. Choose **Public**
5. **Do not** add a README, .gitignore, or license (this project already has them)
6. Click **Create repository**

## 2. Add the remote and push

Replace `YOUR_USERNAME` with your GitHub username (or your org name). If you chose a different repo name, use it instead of `MetagenomeGenerator`.

**HTTPS:**
```bash
cd /home/alex/CursorProjects/MetagenomeGenerator
git remote add origin https://github.com/YOUR_USERNAME/MetagenomeGenerator.git
git push -u origin main
```

**SSH (if you use SSH keys):**
```bash
cd /home/alex/CursorProjects/MetagenomeGenerator
git remote add origin git@github.com:YOUR_USERNAME/MetagenomeGenerator.git
git push -u origin main
```

## 3. Optional: install GitHub CLI for next time

```bash
sudo apt install gh
gh auth login
```

Then you can create and push in one step:  
`gh repo create MetagenomeGenerator --public --source=. --remote=origin --push`
