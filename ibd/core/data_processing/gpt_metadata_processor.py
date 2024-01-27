from openai import OpenAI


class GPTMetadataProcessor():
    def __init__(self):
        self.client = OpenAI()
        
        with open('ibd/core/data_processing/system_prompt.txt', 'r') as f:
            self.system_prompt = f.read()
        

    def process(self, dataset):
        completion = self.client.chat.completions.create(
          model="gpt-4-turbo-preview",
          messages=[
            {"role": "system", "content": self.system_prompt},
            {"role": "user", "content": dataset.to_csv()}
          ],
          temperature=0,
          n=1
        )

        return completion.choices[0].message.content
